# app.py

from flask import Flask, render_template, request, flash
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
app.secret_key = "supersecretkey"  # Replace in production


# -------------------------------------------------
# Helper: Format citation
# -------------------------------------------------
def format_citation(article):
    """Create a human-readable citation string."""
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
# Cache-enabled Search
# -------------------------------------------------
def search_with_cache(source, query, filters):
    """Check cache for results, otherwise fetch from API and save to history."""
    filters_str = json.dumps(filters, sort_keys=True) if isinstance(filters, (dict, list)) else str(filters)

    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT id FROM searches
            WHERE source=? AND query=? AND filters=?
        """, (source, query, filters_str))
        search = cursor.fetchone()

        # ✅ Cache hit
        if search:
            search_id = search[0]
            print(f"[CACHE] ✅ Found results in cache (Search ID={search_id})")

            # Fetch cached articles
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

            # 🔄 Save cached results as a *new* history entry
            save_search(source, query, filters, cached_results)

            return cached_results

    # ❌ Cache miss → Fetch from API
    print("[CACHE] ❌ Not found in cache → Fetching from API")
    results = fetch_from_api(source, query, filters)

    # ✅ Save fresh API results to database
    if results:
        save_search(source, query, filters, results)

    return results


# -------------------------------------------------
# API Search Dispatcher
# -------------------------------------------------
def fetch_from_api(source, query, filters):
    """Route API search to correct client based on source."""
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
@app.route("/", methods=["GET", "POST"])
def home():
    """Home page: Search form + results."""
    results = []
    source = ""
    query = ""
    status_filter = None

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()

        if source == "clinicaltrials":
            status_filter = request.form.get("status", "").strip() or None

        if not source or not query:
            flash("⚠ Please select a source and enter a search query.", "error")
            return render_template("index.html", results=[], source=source, query=query)

        try:
            filters = {"status": status_filter} if source == "clinicaltrials" else {}
            results = search_with_cache(source, query, filters)

            for article in results:
                article["citation"] = format_citation(article)

        except Exception as e:
            flash(f"❌ Error: {str(e)}", "error")

    return render_template(
        "index.html",
        results=results,
        source=source,
        query=query,
        status=status_filter
    )


@app.route("/history")
def history():
    """Display search history."""
    try:
        history_data = get_archive()
    except Exception as e:
        flash(f"⚠ Could not load history: {str(e)}", "error")
        history_data = []

    return render_template("history.html", history=history_data)


@app.route("/history/<int:search_id>")
def history_results(search_id):
    """Display results for a specific history entry."""
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

