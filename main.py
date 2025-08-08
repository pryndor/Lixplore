# main.py

from flask import Flask, render_template, request, flash
from pubmed_client import search_pubmed
from crossref.crossref_client import search_crossref
from springer.springer_client import search_springer
from europepmc.europepmc_client import search_epmc_publications, search_epmc_grants
from clinicaltrials_client import search_clinical_trials
from database.database import get_archive
from database.database import init_db
from doaj_client import search_doaj

# Initialize Flask app
app = Flask(__name__)
app.secret_key = "supersecretkey"  # Required for flash messages


# =====================================================
# Helper function to format citation text
# =====================================================
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


# =====================================================
# Home route - Search interface
# =====================================================
@app.route("/", methods=["GET", "POST"])
def home():
    results = []
    source = ""
    query = ""
    status_filter = None

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()
        status_filter = request.form.get("status", "").strip() or None if source == "clinicaltrials" else None

        if not source or not query:
            flash("⚠ Please select a source and enter a search query.", "error")
            return render_template("index.html", results=[], source=source, query=query)

        try:
            if source == "pubmed":
                results = search_pubmed(query)
            elif source == "crossref":
                results = search_crossref(query)
            elif source == "springer":
                results = search_springer(query)
            elif source == "europepmc_publications":
                results = search_epmc_publications(query)
            elif source == "europepmc_grants":
                grant_data = search_epmc_grants(query)
                results = grant_data.get("grants", []) if isinstance(grant_data, dict) else []
            elif source == "clinicaltrials":
                results = search_clinical_trials(query, status_filter)
            elif source == "doaj":
                results = search_doaj(query)
            else:
                flash("⚠ Unknown source selected.", "error")

            for article in results:
                article["citation"] = format_citation(article)

        except Exception as e:
            flash(f"❌ Error while searching: {str(e)}", "error")
            results = []

    return render_template("index.html", results=results, source=source, query=query, status=status_filter)


# =====================================================
# History route - Previously stored searches
# =====================================================
@app.route("/history")
def history():
    try:
        history_data = get_archive()
    except Exception as e:
        flash(f"⚠ Could not load history: {str(e)}", "error")
        history_data = []

    return render_template("history.html", history=history_data)


# =====================================================
# App entry point
# =====================================================
if __name__ == "__main__":
    print("🔧 Initializing blog database...")
    init_db()
    print("✅ Blog database initialized.")
    app.run(host="0.0.0.0", port=5000, debug=True)

