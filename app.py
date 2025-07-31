#app.py

from flask import Flask, render_template, request, flash
from pubmed_client import search_pubmed
from crossref.crossref_client import search_crossref
from springer.springer_client import search_springer
from europepmc.europepmc_client import search_epmc_publications, search_epmc_grants
from database.database import get_archive

app = Flask(__name__)
app.secret_key = "supersecretkey"  # Required for flash messages

@app.route("/", methods=["GET", "POST"])
def home():
    results = []

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()

        # Validate inputs
        if not source or not query:
            flash("⚠ Please select a source and enter a search query.", "error")
            return render_template("index.html", results=[])

        try:
            # Route search based on source
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
                # Europe PMC grants may need custom display
                results = grant_data.get("grants", []) if isinstance(grant_data, dict) else []

            else:
                flash("⚠ Unknown source selected.", "error")

        except Exception as e:
            flash(f"❌ An error occurred while searching: {str(e)}", "error")
            results = []

    return render_template("index.html", results=results)


@app.route("/history")
def history():
    try:
        history_data = get_archive()
    except Exception as e:
        flash(f"⚠ Could not load history: {str(e)}", "error")
        history_data = []

    return render_template("history.html", history=history_data)


if __name__ == "__main__":
    app.run(debug=True)

