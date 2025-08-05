from flask import Flask, render_template, request, flash, jsonify
from pubmed_client import search_pubmed
from crossref.crossref_client import search_crossref
from springer.springer_client import search_springer
from europepmc.europepmc_client import search_epmc_publications, search_epmc_grants
from clinicaltrials_client import search_clinical_trials
from doaj_client import search_doaj
from database.database import (
    init_db,
    get_archive,
    get_connection,
    save_search,
    get_results_by_search_id
)
import json

# -------------------------------------------------
# Flask App Configuration
# -------------------------------------------------
app = Flask(__name__)
app.secret_key = "supersecretkey"  # 🔴 Replace with a secure key in production


# -------------------------------------------------
# Helper: Format citation
# -------------------------------------------------
def format_citation(article):
    authors = ", ".join(article.get("authors", [])) if article.get("authors") else "Unknown Author"
    year = article.get("year", "n.d.")
    title = article.get("title", "No title available")
    doi = article.get("doi")
    url = article.get("url")

    citation = f"{authors} ({year}). {title}."
    if doi:
        citation += f" https://doi.org/{doi}"
    elif url:
        citation += f" {url}"
    return citation


# -------------------------------------------------
# Search Dispatcher with Cache
# -------------------------------------------------
def search_with_cache(source, query, filters):
    filters_str = json.dumps(filters, sort_keys=True) if isinstance(filters, (dict, list)) else str(filters)

    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT id FROM searches
            WHERE source=? AND query=? AND filters=?
        """, (source, query, filters_str))
        search = cursor.fetchone()

        if search:
            search_id = search[0]
            rows = get_results_by_search_id(search_id)
            cached_results = [
                {
                    "title": row[0],
                    "authors": row[1].split(", ") if row[1] else [],
                    "doi": row[2],
                    "pmid": row[3],
                    "url": row[4],
                    "abstract": row[5],
                    "affiliations": row[6].split(", ") if row[6] else [],
                    "oa_pdf_path": row[7]
                }
                for row in rows
            ]
            save_search(source, query, filters, cached_results)
            return cached_results

    # If not cached, fetch from API
    results = fetch_from_api(source, query, filters)
    if results:
        save_search(source, query, filters, results)
    return results


def fetch_from_api(source, query, filters):
    try:
        if source == "pubmed":
            return search_pubmed(query)
        elif source == "crossref":
            return search_crossref(query)
        elif source == "springer":
            return search_springer(query)
        elif source == "europepmc_publications":
            return search_epmc_publications(query)
        elif source == "europepmc_grants":
            grant_data = search_epmc_grants(query)
            return grant_data.get("grants", []) if isinstance(grant_data, dict) else []
        elif source == "clinicaltrials":
            return search_clinical_trials(query, filters.get("status"))
        elif source == "doaj":
            return search_doaj(query)
    except Exception as e:
        print(f"[API] ⚠ Error fetching from {source}: {e}")
    return []


# -------------------------------------------------
# Routes
# -------------------------------------------------

# Intro Page
@app.route("/")
def intro():
    return render_template("intro.html")


# Search Page
@app.route("/search", methods=["GET", "POST"])
def search_page():
    results, source, query, status_filter = [], "", "", None

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()
        status_filter = request.form.get("status", "").strip() or None if source == "clinicaltrials" else None

        if not source or not query:
            flash("⚠ Please select a source and enter a query.", "error")
            return render_template("index.html", results=[], source=source, query=query)

        try:
            filters = {"status": status_filter} if status_filter else {}
            results = search_with_cache(source, query, filters)
            for article in results:
                article["citation"] = format_citation(article)
        except Exception as e:
            flash(f"❌ Error: {str(e)}", "error")

    return render_template("index.html", results=results, source=source, query=query, status=status_filter)


# Zotero Search (GET Parameters)
@app.route("/zotero_search", methods=["GET"])
def zotero_search():
    title = request.args.get("title", "")
    author = request.args.get("author", "")
    year = request.args.get("year", "")
    doi = request.args.get("doi", "")

    search_term = doi or title or author

    if not search_term:
        return jsonify({"error": "No valid search parameters"}), 400

    results = search_with_cache("crossref", search_term, {})
    for article in results:
        article["citation"] = format_citation(article)

    return render_template("index.html", results=results, source="crossref", query=search_term)


# History Page
@app.route("/history")
def history():
    try:
        history_data = get_archive()
    except Exception as e:
        flash(f"⚠ Could not load history: {str(e)}", "error")
        history_data = []
    return render_template("history.html", history=history_data)


# History Results Page
@app.route("/history/<int:search_id>")
def history_results(search_id):
    try:
        rows = get_results_by_search_id(search_id)
        articles = [
            {
                "title": row[0],
                "authors": row[1].split(", ") if row[1] else [],
                "doi": row[2],
                "pmid": row[3],
                "url": row[4],
                "abstract": row[5],
                "affiliations": row[6].split(", ") if row[6] else [],
                "oa_pdf_path": row[7],
                "citation": format_citation({
                    "authors": row[1].split(", ") if row[1] else [],
                    "year": "n.d.",
                    "title": row[0],
                    "doi": row[2],
                    "url": row[4]
                })
            }
            for row in rows
        ]
    except Exception as e:
        flash(f"⚠ Could not load search results: {str(e)}", "error")
        articles = []
    return render_template("results.html", results=articles, source="History Search")


# -------------------------------------------------
# Main
# -------------------------------------------------
if __name__ == "__main__":
    init_db()
    app.run(host="0.0.0.0", port=5000, debug=True)

