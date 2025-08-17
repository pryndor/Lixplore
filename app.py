# app.py

from flask import Flask, render_template, request, flash, jsonify, redirect, url_for
import markdown
import sqlite3
import os
import json
from datetime import datetime
import re
import yaml
import smtplib

from cache import get_cached_results, save_results_to_cache
from pubmed_client import search_pubmed
from crossref.crossref_client import search_crossref
from springer.springer_client import search_springer
from europepmc.europepmc_client import search_epmc_publications, search_epmc_grants
from clinicaltrials_client import search_clinical_trials
from doaj_client import search_doaj
from flask_mail import Mail, Message

from database.database import (
    init_db,
    get_archive,
    DB_PATH,
    get_connection,
    save_search,
    get_results_by_search_id,
    get_posts_by_tag
)

app = Flask(__name__)
app.secret_key = "supersecretkey"  # Change in production
DATABASE = os.path.join(os.path.dirname(__file__), "database", "lixplore.db")


def get_db_connection():
    conn = sqlite3.connect(DATABASE)
    conn.row_factory = sqlite3.Row  # Allows dict-like access to rows
    return conn


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

def save_search_and_articles(source, query, filters, results):
    filters_str = json.dumps(filters or {}, sort_keys=True)

    with get_connection() as conn:
        cursor = conn.cursor()

        # Insert or get search_id
        cursor.execute("""
            SELECT id FROM searches WHERE source=? AND query=? AND filters=?
        """, (source, query, filters_str))
        row = cursor.fetchone()
        if row:
            search_id = row[0]
        else:
            cursor.execute("""
                INSERT INTO searches (source, query, filters) VALUES (?, ?, ?)
            """, (source, query, filters_str))
            search_id = cursor.lastrowid

        # Clear old articles for this search if any
        cursor.execute("DELETE FROM articles WHERE search_id=?", (search_id,))

        # Insert articles individually
        for article in results:
            authors_str = ", ".join(article.get("authors", [])) if article.get("authors") else None
            affiliations_str = ", ".join(article.get("affiliations", [])) if article.get("affiliations") else None

            cursor.execute("""
                INSERT INTO articles (
                    search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                search_id,
                article.get("title"),
                authors_str,
                article.get("doi"),
                article.get("pmid"),
                article.get("url"),
                article.get("abstract"),
                affiliations_str,
                article.get("oa_pdf_path"),
            ))

        # Optionally update cache table or ignore it if articles hold the data
        # Commit transaction
        conn.commit()

def get_results_by_search_id(search_id):
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
            FROM articles
            WHERE search_id=?
        """, (search_id,))
        return cursor.fetchall()

def parse_front_matter(content):
    """
    Splits the Markdown front matter (YAML) from the body.
    Returns: (metadata_dict, body_text)
    """
    match = re.match(r'^---\n(.*?)\n---\n(.*)', content, re.S)
    if match:
        front_matter = match.group(1)
        body = match.group(2)
        metadata = yaml.safe_load(front_matter)
        return metadata, body
    else:
        return {}, content


def load_markdown_posts():
    conn = get_db_connection()
    cursor = conn.cursor()

    posts_dir = os.path.join(os.path.dirname(__file__), "posts")  # adjust if needed

    for filename in os.listdir(posts_dir):
        if filename.endswith(".md"):
            filepath = os.path.join(posts_dir, filename)
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

            # Parse front matter
            metadata, body = parse_front_matter(content)

            # Convert the markdown body to HTML
            body_html = markdown.markdown(body)
            title = metadata.get('title', 'Untitled')
            author = metadata.get('author', 'Bala PR')  # <-- your default
            published = metadata.get('published', datetime.now().strftime("%Y-%m-%d"))
            tags = metadata.get('tags', '')
            slug = metadata.get('slug', filename[:-3])

            # Check if already in DB
            cursor.execute("SELECT COUNT(*) FROM articles WHERE slug = ?", (slug,))
            if cursor.fetchone()[0] == 0:
                cursor.execute(
                    "INSERT INTO articles (title, author, published, tags, slug, content) VALUES (?, ?, ?, ?, ?, ?)",
                    (title, author, published, tags, slug, body_html)
                )

    conn.commit()
    conn.close()




def get_cached_results(source, query, filters):
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
            cursor.execute("""
                SELECT title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                FROM results WHERE search_id=?
            """, (search_id,))
            rows = cursor.fetchall()
            results = [
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
            return results, search_id
    return None, None


import json

def safe_str(value):
    """Convert any value to a string for DB storage. Dicts/lists -> JSON string, None -> None"""
    if value is None:
        return None
    if isinstance(value, str):
        return value
    try:
        return json.dumps(value)
    except Exception:
        return str(value)


def save_results_to_cache(source, query, filters, results):
    filters_str = json.dumps(filters or {}, sort_keys=True)
    with get_connection() as conn:
        cursor = conn.cursor()

        # Insert search entry
        cursor.execute("""
            INSERT INTO searches (source, query, filters) VALUES (?, ?, ?)
        """, (source, query, filters_str))
        search_id = cursor.lastrowid

        # Insert each article
        for article in results:
            authors = safe_str(article.get("authors"))
            affiliations = safe_str(article.get("affiliations"))
            abstract = safe_str(article.get("abstract"))
            oa_pdf_url = article.get("oa_pdf_path")
            if not isinstance(oa_pdf_url, str):
                oa_pdf_url = safe_str(oa_pdf_url)

            cursor.execute("""
                INSERT INTO results (
                    search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                search_id,
                safe_str(article.get("title")),
                authors,
                safe_str(article.get("doi")),
                safe_str(article.get("pmid")),
                safe_str(article.get("url")),
                abstract,
                affiliations,
                oa_pdf_url
            ))

        conn.commit()



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

#--------------------------
# Mail config----
#------------------------

# Configure your email server (example uses Gmail SMTP)
app.config.update(
    MAIL_SERVER='smtppro.zoho.in',
    MAIL_PORT=587,
    MAIL_USE_TLS=True,
    MAIL_USE_SSL=False,
    MAIL_USERNAME='info@drugvigil.com',        # your sending email
    MAIL_PASSWORD='pyJECPvDkLs6',  # password or app password
    MAIL_DEFAULT_SENDER=('Lixlpore Alerts', 'info@drugvigil.com')
)

mail = Mail(app)


# -------------------------------------------------
# Routes
# -------------------------------------------------


@app.route("/contact", methods=["GET", "POST"])
def contact():
    if request.method == "POST":
        name = request.form.get("name")
        email = request.form.get("email")
        subject = request.form.get("subject")
        message = request.form.get("message")

        if not (name and email and subject and message):
            flash("Please fill in all the fields.", "error")
            return redirect(url_for("contact"))

        # Compose email content
        email_body = f"""
        New message from Contact Form:

        Name: {name}
        Email: {email}
        Subject: {subject}

        Message:
        {message}
        """

        try:
            msg = Message(
                subject=f"Contact Form: {subject}",
                sender=(name, email),  # Name and email of sender
                recipients=["info@drugvigil.com"],  # Replace with your receiving email
                body=email_body
            )
            mail.send(msg)
            flash("Message sent successfully!", "success")
        except Exception as e:
            flash(f"Failed to send message: {e}", "error")

        return redirect(url_for("contact"))

    # If GET request, show the contact form
    return render_template("contact.html")

    

@app.route("/blogs/tag/<tag>")
def blogs_by_tag(tag):
    conn = get_db_connection()
    posts = conn.execute("""
        SELECT title, author, created_at, tags, content, slug
        FROM blog_posts
        WHERE tags LIKE ?
        ORDER BY created_at DESC
    """, ('%' + tag + '%',)).fetchall()
    conn.close()

    # Convert sqlite3.Row to dict (optional, depends on how you pass to template)
    posts = [dict(post) for post in posts]

    return render_template("blog.html", posts=posts)

@app.route("/")
def intro():
    return render_template("intro.html")

@app.route("/blog")
def blog():
    conn = get_db_connection()
    posts = conn.execute("""
        SELECT title, author, created_at, tags, content, slug
        FROM blog_posts
        ORDER BY created_at DESC
    """).fetchall()
    conn.close()

    posts_list = []
    for post in posts:
        html_content = markdown.markdown(post["content"])  # Convert markdown to HTML here
        posts_list.append({
            "title": post["title"],
            "author": post["author"] or "Admin",
            "created_at": post["created_at"],
            "tags": post["tags"] or "",
            "content": html_content,
            "slug": post["slug"]
        })
    
    return render_template("blog.html", posts=posts_list)


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
             # Test PMC PDF URL
            print("PMC PDF URL:", article.get("pmc_pdf_url"))


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

