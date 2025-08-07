# app.py
from flask import Flask, render_template, request, flash, jsonify
import markdown
import os
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

app = Flask(__name__)
app.secret_key = "supersecretkey"  # ⚠ Change in production


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
# Helper: Search with cache
# -------------------------------------------------
def search_with_cache(source, query, filters):
    filters_str = json.dumps(filters or {}, sort_keys=True)

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
            return [
                {
                    "title": r[0],
                    "authors": r[1].split(", ") if r[1] else [],
                    "doi": r[2],
                    "pmid": r[3],
                    "url": r[4],
                    "abstract": r[5],
                    "affiliations": r[6].split(", ") if r[6] else [],
                    "oa_pdf_path": r[7]
                }
                for r in rows
            ]

    # Not cached → fetch fresh
    results = fetch_from_api(source, query, filters)
    if results:
        save_search(source, query, filters, results)
    return results


# -------------------------------------------------
# Helper: API fetcher
# -------------------------------------------------
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
@app.route("/")
def intro():
    return render_template("intro.html")

@app.route("/blog")
def blog():
    blog_dir = "blog_posts"
    posts = {}

    for filename in sorted(os.listdir(blog_dir), reverse=True):
        if filename.endswith(".md"):
            filepath = os.path.join(blog_dir, filename)
            with open(filepath, "r", encoding="utf-8") as f:
                md_content = f.read()
                html_content = markdown.markdown(md_content)
                title = filename.replace(".md", "").replace("-", " ").title()
                posts[title] = html_content

    return render_template("blog.html", posts=posts)


@app.route("/release-notes")
def show_release_notes():
    base_dir = "release_notes"
    releases = {}

    if not os.path.exists(base_dir):
        return "Release notes directory not found", 404

    for folder in sorted(os.listdir(base_dir), reverse=True):
        folder_path = os.path.join(base_dir, folder)

        # Ensure only directories are considered
        if not os.path.isdir(folder_path):
            continue

        note_file = os.path.join(folder_path, "notes.md")

        if os.path.exists(note_file):
            with open(note_file, "r", encoding="utf-8") as f:
                md_content = f.read()
                html_content = markdown.markdown(md_content)
                releases[folder] = html_content

    print("Releases loaded:", releases.keys())


    return render_template("release_notes.html", releases=releases)


@app.route("/load_more", methods=["GET"])
def load_more():
    try:
        source = request.args.get("source")
        query = request.args.get("query", "").strip()
        start = int(request.args.get("start", 0))
        max_results = int(request.args.get("max", 20))

        if not source or not query:
            return jsonify({"error": "Missing source or query"}), 400

        # Call unified fetcher (handles all sources)
        results = fetch_from_api(source, query, {"start": start, "max": max_results}) or []

        # Add formatted citations
        for article in results:
            article["citation"] = format_citation(article)

        return jsonify({"results": results, "next_start": start + max_results})

    except Exception as e:
        return jsonify({"error": str(e)}), 500



@app.route("/privacy")
def privacy():
    return render_template("pandt.html")


@app.route("/support")
def support():
    return render_template("support.html")


@app.route("/search", methods=["GET", "POST"])
def search_page():
    results, source, query, status_filter = [], "", "", None

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()
        status_filter = request.form.get("status", "").strip() if source == "clinicaltrials" else None

        if not source or not query:
            flash("⚠ Please select a source and enter a query.", "error")
            return render_template("index.html", results=[], source=source, query=query)

        filters = {"status": status_filter} if status_filter else {}
        results = search_with_cache(source, query, filters)

        for article in results:
            article["citation"] = format_citation(article)

    return render_template("index.html", results=results, source=source, query=query, status=status_filter)


@app.route("/history")
def history():
    history_data = get_archive()
    return render_template("history.html", history=history_data)


@app.route("/donations")
def donations():
    return render_template("donations.html")


@app.route("/blogs")
def blogs():
    return render_template("blogs.html")


@app.route("/contact")
def contact():
    return render_template("contact.html")


@app.route("/history/<int:search_id>")
def history_results(search_id):
    rows = get_results_by_search_id(search_id)
    articles = [
        {
            "title": r[0],
            "authors": r[1].split(", ") if r[1] else [],
            "doi": r[2],
            "pmid": r[3],
            "url": r[4],
            "abstract": r[5],
            "affiliations": r[6].split(", ") if r[6] else [],
            "oa_pdf_path": r[7],
            "citation": format_citation({
                "authors": r[1].split(", ") if r[1] else [],
                "year": "n.d.",
                "title": r[0],
                "doi": r[2],
                "url": r[4]
            })
        }
        for r in rows
    ]
    return render_template("results.html", results=articles, source="History Search")


# -------------------------------------------------
# Main Entry
# -------------------------------------------------
if __name__ == "__main__":
    init_db()
    app.run(host="0.0.0.0", port=5000, debug=True)

