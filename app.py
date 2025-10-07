# app.py

from flask import Flask, render_template, request, flash, jsonify, redirect, url_for, session
import markdown
import sqlite3
import os
import json
from datetime import datetime
import re
import yaml
import smtplib
import requests
from markupsafe import escape

from cache import get_cached_results, save_results_to_cache
from api_clients.pubmed.pubmed_client import search_pubmed
from api_clients.crossref.crossref_client import search_crossref
from api_clients.springer.springer_client import search_springer
from api_clients.europepmc.europepmc_client import search_epmc_publications, search_epmc_grants
from api_clients.clinicaltrial.clinicaltrials_client import search_clinical_trials
from api_clients.doaj.doaj_client import search_doaj
from flask_mail import Mail, Message



from database.database import (
    init_db,
    get_archive,
    get_search_history_by_type,
    get_matched_fields,
    DB_PATH,
    get_connection,
    save_search,
    get_results_by_search_id,
    get_posts_by_tag,
    # Screening functions
   # create_screening_project,
   # get_screening_projects,
   # get_screening_project,
   # add_reviewer_to_project,
   # get_project_reviewers,
   # add_articles_to_screening_project,
   # get_screening_articles,
   # save_screening_decision,
   # get_screening_statistics
)

app = Flask(__name__)
app.secret_key = "supersecretkey"  # Change in production
DATABASE = os.path.join(os.path.dirname(__file__), "database", "lixplore.db")


def get_db_connection():
    conn = sqlite3.connect(DATABASE)
    conn.row_factory = sqlite3.Row  # Allows dict-like access to rows
    return conn


# -------------------------------------------------
# Helper: API fetcher
# -------------------------------------------------

def fetch_from_api(source, query, filters):
    try:
        if source == "pubmed":
            results_dict = search_pubmed(query)
            print("Query received:", query)
            results = results_dict.get("results", [])

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
            results = search_clinical_trials(query, filters.get("status"))
        elif source == "doaj":
            results = search_doaj(query)
        else:
            results = []

        # Normalize OA PDF paths for all articles
        for i, article in enumerate(results):
            results[i] = normalize_oa_pdf(article)

        return results

    except Exception as e:
        print(f"[API] ⚠ Error fetching from {source}: {e}")
        return []



def get_search_results():
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute("SELECT * FROM articles ORDER BY id DESC LIMIT 20")
    rows = cur.fetchall()
    results = [dict(row) for row in rows]
    conn.close()
    return results



# -----------------------------
# Helper: Normalize OA PDF
# -----------------------------



def normalize_oa_pdf(article):
    oa_pdf = article.get("oa_pdf_path")
    if isinstance(oa_pdf, dict):
        article["oa_pdf_path"] = oa_pdf.get("value")
    elif isinstance(oa_pdf, list):
        article["oa_pdf_path"] = oa_pdf[0] if oa_pdf else None
    return article

#------------------------
# Tagging function under implementation: functions under testing
#------------------------


# Add tags to article
def add_tags_to_article(article_id, tag_names):
    conn = get_db_connection()
    cur = conn.cursor()
    for tag_name in tag_names:
        cur.execute("INSERT OR IGNORE INTO article_tags (name) VALUES (?)", (tag_name.strip(),))
        cur.execute("SELECT id FROM article_tags WHERE name = ?", (tag_name.strip(),))
        tag_id = cur.fetchone()["id"]
        cur.execute("INSERT OR IGNORE INTO article_tag_map (article_id, tag_id) VALUES (?, ?)", (article_id, tag_id))
    conn.commit()
    conn.close()


# Get tags for article
def get_tags_for_article(article_id):
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute("""
        SELECT t.name
        FROM article_tags t
        JOIN article_tag_map m ON t.id = m.tag_id
        WHERE m.article_id = ?
    """, (article_id,))
    tags = [row["name"] for row in cur.fetchall()]
    conn.close()
    return tags


# Get articles by tag
def get_articles_by_tag(tag_name):
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute("""
        SELECT a.*
        FROM articles a
        JOIN article_tag_map m ON a.id = m.article_id
        JOIN article_tags t ON m.tag_id = t.id
        WHERE t.name = ?
    """, (tag_name,))
    articles = [dict(row) for row in cur.fetchall()]  # always return dicts
    conn.close()
    return articles


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



# -----------------------------
# Search with cache
# -----------------------------
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

    # Normalize OA PDF before saving
    for i, article in enumerate(results):
        results[i] = normalize_oa_pdf(article)

    if results:
        save_search_and_articles(source, query, filters, results)

    return results

# -----------------------------
# Save search and articles
# -----------------------------
def save_search_and_articles(source, query, filters, results):
    filters_str = json.dumps(filters or {}, sort_keys=True)
    with get_connection() as conn:
        cursor = conn.cursor()

        # Insert or get search_id
        cursor.execute("SELECT id FROM searches WHERE source=? AND query=? AND filters=?",
                       (source, query, filters_str))
        row = cursor.fetchone()
        if row:
            search_id = row[0]
        else:
            cursor.execute("INSERT INTO searches (source, query, filters) VALUES (?, ?, ?)",
                           (source, query, filters_str))
            search_id = cursor.lastrowid

        # Clear old articles
        cursor.execute("DELETE FROM articles WHERE search_id=?", (search_id,))

        for article in results:
            article = normalize_oa_pdf(article)  # Ensure OA PDF is string
            authors_str = ", ".join(article.get("authors", [])) if article.get("authors") else None
            affiliations_str = ", ".join(article.get("affiliations", [])) if article.get("affiliations") else None

            cursor.execute("""
                INSERT INTO articles (search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                search_id,
                safe_str(article.get("title")),
                safe_str(authors_str),
                safe_str(article.get("doi")),
                safe_str(article.get("pmid")),
                safe_str(article.get("url")),
                safe_str(article.get("abstract")),
                safe_str(affiliations_str),
                safe_str(article.get("oa_pdf_path"))
            ))
        conn.commit()
'''
def get_results_by_search_id(search_id):
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
            FROM articles
            WHERE search_id=?
        """, (search_id,))
        columns = [desc[0] for desc in cursor.description]
        rows = cursor.fetchall()

        # Convert rows -> dicts (no normalization here)
        results = [dict(zip(columns, row)) for row in rows]

        return results
'''



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




# -----------------------------
# Helper: Safe string conversion
# -----------------------------
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


# -----------------------------
# Save results to cache (optional)
# -----------------------------
def save_results_to_cache(source, query, filters, results):
    filters_str = json.dumps(filters or {}, sort_keys=True)
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("INSERT INTO searches (source, query, filters) VALUES (?, ?, ?)",
                           (source, query, filters_str))
            search_id = cursor.lastrowid

            for article in results:
                article = normalize_oa_pdf(article)  # Normalize before saving
                authors_str = ", ".join(article.get("authors", [])) if article.get("authors") else None
                affiliations_str = ", ".join(article.get("affiliations", [])) if article.get("affiliations") else None

                cursor.execute("""
                    INSERT INTO articles (
                        search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    search_id,
                    safe_str(article.get("title")),
                    safe_str(authors_str),
                    safe_str(article.get("doi")),
                    safe_str(article.get("pmid")),
                    safe_str(article.get("url")),
                    safe_str(article.get("abstract")),
                    safe_str(affiliations_str),
                    safe_str(article.get("oa_pdf_path"))
                ))
                print("DEBUG PDF:", article.get("oa_pdf_path"))  # 👈 add this
            conn.commit()
    except Exception as e:
        print(f"[DB] ⚠ Error saving search: {e}")


# -------------------------------------------------
# Helper: API fetcher
# -------------------------------------------------


"""
def fetch_from_api(source, query, filters):
    try:
        if source == "pubmed":
            results = search_pubmed(query)
            print("Query received:", query)

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
            results = search_clinical_trials(query, filters.get("status"))
        elif source == "doaj":
            results = search_doaj(query)
        else:
            results = []

        # Normalize OA PDF paths for all articles
        for i, article in enumerate(results):
            results[i] = normalize_oa_pdf(article)

        return results

    except Exception as e:
        print(f"[API] ⚠ Error fetching from {source}: {e}")
        return []


"""
# -----------------------------
# Helper: Normalize OA PDF
# -----------------------------
def normalize_oa_pdf(article):
    oa_pdf = article.get("oa_pdf_path")
    if isinstance(oa_pdf, dict):
        article["oa_pdf_path"] = oa_pdf.get("value")
    elif isinstance(oa_pdf, list):
        article["oa_pdf_path"] = oa_pdf[0] if oa_pdf else None
    return article

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

@app.route('/home')
def home():
    return render_template('index.html')  # Or your main search page template

@app.template_filter('escapejs')
def escapejs_filter(s):
    if s is None:
        return ''
    return json.dumps(s)


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

        # Compose email content (dedent to avoid spacing issues)
        email_body = textwrap.dedent(f"""\
            New message from Contact Form:

            Name: {name}
            Email: {email}
            Subject: {subject}

            Message:
            {message}
        """)

        try:
            msg = Message(
                subject=f"Contact Form: {subject}",
                sender=(name, email),   # sender tuple: (name, email)
                recipients=["info@drugvigil.com"],  # change to your receiving email
                body=email_body
            )
            mail.send(msg)
            flash("Message sent successfully!", "success")
        except Exception as e:
            flash(f"Failed to send message: {e}", "error")

        return redirect(url_for("contact"))

    # GET request → render contact form
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


# Load more functionality removed - now using pagination


@app.route("/privacy")
def privacy():
    return render_template("pandt.html")


@app.route("/support")
def support():
    return render_template("support.html")

"""
@app.route("/search", methods=["GET", "POST"])
def search_page():
    results, source, query, status_filter = [], "", "", None

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()
        status_filter = request.form.get("status", "").strip() if source == "clinicaltrials" else None

    elif request.method == "GET" and request.args.get('source') and request.args.get('query'):
        source = request.args.get("source")
        query = request.args.get("query").strip()
        status_filter = request.args.get("status")

    # Handle search only if source + query are provided
    if source and query:
        filters = {"status": status_filter} if status_filter else {}
        results = search_with_cache(source, query, filters)

        # Sort results
        results.sort(key=lambda x: x.get("title", "").lower())

        # Post-process results (shared)
        for article in results:
            pdf_url = article.get("pmc_pdf_url") or article.get("oa_pdf_path")
            article["oa_pdf_path"] = pdf_url
            article["citation"] = format_citation(article)
    else:
        if request.method == "POST":  # only warn on invalid POST
            flash("⚠ Please select a source and enter a query.", "error")

    return render_template(
        "index.html",
        results=results,
        source=source,
        query=query,
        status=status_filter
    )


"""
@app.route("/search", methods=["GET", "POST"])
def search_page():
    results, source, query, status_filter = [], "", "", None

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()
        status_filter = request.form.get("status", "").strip() if source == "clinicaltrials" else None

        if not source or not query:
            flash("⚠ Please select a source and enter a query.", "error")
            return render_template("index.html", results=[], source=source, query=query, status=status_filter)

        filters = {"status": status_filter} if status_filter else {}
        results = search_with_cache(source, query, filters)

        # Sort results alphabetically by title
        results.sort(key=lambda x: x.get("title", "").lower())

        # Process results
        for article in results:
            pdf_url = article.get("pmc_pdf_url") or article.get("oa_pdf_path")
            article["oa_pdf_path"] = pdf_url
            article["citation"] = format_citation(article)

    elif request.method == "GET" and request.args.get('source') and request.args.get('query'):
        source = request.args.get("source")
        query = request.args.get("query")
        status_filter = request.args.get("status")

        filters = {"status": status_filter} if status_filter else {}
        results = search_with_cache(source, query, filters)
        results.sort(key=lambda x: x.get("title", "").lower())
        for article in results:
            pdf_url = article.get("pmc_pdf_url") or article.get("oa_pdf_path")
            article["oa_pdf_path"] = pdf_url
            article["citation"] = format_citation(article)

    return render_template(
        "index.html",
        results=results,
        source=source,
        query=query,
        status=status_filter
    )



# -----------------------------
# Helper function: get tags for an article
# -----------------------------
def get_tags_for_article(article_id):
    conn = sqlite3.connect("database/lixplore.db")
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute("""
        SELECT t.name
        FROM article_tags t
        JOIN article_tag_map m ON t.id = m.tag_id
        WHERE m.article_id = ?
    """, (article_id,))
    tags = [row["name"] for row in cur.fetchall()]
    conn.close()
    return tags

# -----------------------------
# Helper function: add tags to an article
# -----------------------------
def add_tags_to_article(article_id, tags):
    conn = sqlite3.connect("database/lixplore.db")
    cur = conn.cursor()
    for tag_name in tags:
        tag_name = tag_name.strip()
        if not tag_name:
            continue
        # Check if tag already exists in article_tags
        cur.execute("SELECT id FROM article_tags WHERE name = ?", (tag_name,))
        tag_row = cur.fetchone()
        if tag_row:
            tag_id = tag_row[0]
        else:
            # Insert new tag
            cur.execute("INSERT INTO article_tags (name) VALUES (?)", (tag_name,))
            tag_id = cur.lastrowid
        # Map tag to article
        cur.execute("SELECT * FROM article_tag_map WHERE article_id = ? AND tag_id = ?", (article_id, tag_id))
        if not cur.fetchone():
            cur.execute("INSERT INTO article_tag_map (article_id, tag_id) VALUES (?, ?)", (article_id, tag_id))
    conn.commit()
    conn.close()

# -----------------------------
# Route: Search results page
# -----------------------------
"""
@app.route("/search_results")
def search_results():
    results = get_search_results()  # replace with your DB query
    # Ensure each article has id, title, abstract, authors, etc.
    for article in results:
        if 'id' not in article:
            raise ValueError(f"Article missing ID: {article}")
        article['tags'] = get_tags_for_article(article['id'])
    return render_template("results.html", results=results, source="Search")

"""

@app.route("/search_results")
def search_results():
    results = get_search_results()  # your DB query
    print("Results returned:", len(results))
    for article in results:
        print("Article ID:", article.get("id"), "Title:", article.get("title"))
        article['tags'] = get_tags_for_article(article['id'])
    return render_template("results.html", results=results, source="Search")


# -----------------------------
# Route: Add tag to an article
# -----------------------------
@app.route("/add_tags/<int:article_id>", methods=["POST"])
def add_tags(article_id):
    tags_text = request.form.get("tags", "")
    tags = [t.strip() for t in tags_text.split(",") if t.strip()]
    if tags:
        add_tags_to_article(article_id, tags)
    # After adding, redirect back to search results
    return redirect(url_for("search_results"))


@app.route("/articles_with_tags")
def articles_with_tags():
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute("""
        SELECT a.id, a.title, GROUP_CONCAT(t.name, ', ') AS tags
        FROM articles a
        LEFT JOIN article_tag_map m ON a.id = m.article_id
        LEFT JOIN article_tags t ON m.tag_id = t.id
        GROUP BY a.id
        ORDER BY a.id DESC
        LIMIT 20;
    """)
    articles = [dict(row) for row in cur.fetchall()]
    conn.close()
    return render_template("articles_with_tags.html", articles=articles)

#----------------------


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

''' # This is commented out that believed to be unecessary implementaion of pagination for history section, instead implemented filters

@app.route("/search-history")
def search_history():
    """
    Display search results from history without pagination.
    """
    # Retrieve all history from session
    all_history = session.get('search_history', [])

    # Query parameters
    search_type = request.args.get('type', 'local')  # 'local' or 'global'
    search_query = request.args.get('q', '').strip()
    source_filter = request.args.get('source', '')

    # Filter history
    if search_query or source_filter:
        filtered_history = get_search_history_by_type(
            search_type=search_type,
            search_query=search_query,
            source_filter=source_filter
        )
    else:
        filtered_history = get_archive()  # fallback to full archive

    # Sort alphabetically by query
    filtered_history.sort(key=lambda x: x.get('query', '').lower())

    # Annotate each item for match highlights
    for item in filtered_history:
        item['display_type'] = search_type
        item['is_local_match'] = search_query.lower() in item.get('query', '').lower() if search_query else False
        item['is_global_match'] = (
            search_query.lower() in item.get('query', '').lower() or
            search_query.lower() in item.get('source', '').lower() or
            search_query.lower() in str(item.get('timestamp', '')).lower()
        ) if search_query else False

    # Available sources for filtering
    available_sources = sorted({h.get('source', '') for h in all_history if h.get('source')})

    return render_template(
        "history.html",
        history=filtered_history,  # pass all filtered results
        search_query=search_query,
        search_type=search_type,
        source_filter=source_filter,
        available_sources=available_sources
    )

'''

@app.route("/search-history")
def search_history():
    history_data = get_archive()
    return render_template("history.html", history=history_data)

"""
@app.route("/donations")
def donations():
    return render_template("donations.html")


@app.route("/blogs")
def blogs():
    return render_template("blogs.html")

@app.route('/home') #This section updates for database error
def home_page():
    return render_template(
        'index.html',# Or your main search page template
        page=1,
        results=[],
        total=0,
        query=""
    )

"""
#-------------------------
# approute for pagination 
#---------------------------

@app.route('/search/pubmed')
def pubmed_search():
    print("PubMed search route hit")
    query = request.args.get('q', '')
    page = int(request.args.get('page', 1))
    per_page = int(request.args.get('per_page', 20))

    print("Page:", page)
    print("Per Page:", per_page)
    print("Query received:", query)
    return {
        'results': id_list,
        'page': page,
        'per_page': per_page,
        'total': total_count
    }

    # Call the actual search function (you need to ensure this function exists and accepts these params)
    result = search_pubmed(query, page, per_page)

    print("Returning result:", result)

    return jsonify(result)

"""
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

"""

# Pagination and search routes removed as requested


# -------------------------------------------------
# Main Entry
# -------------------------------------------------
if __name__ == "__main__":
    init_db()
    app.run(host="0.0.0.0", port=5000, debug=True)



