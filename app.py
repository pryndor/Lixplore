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
import requests

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
    create_screening_project,
    get_screening_projects,
    get_screening_project,
    add_reviewer_to_project,
    get_project_reviewers,
    add_articles_to_screening_project,
    get_screening_articles,
    save_screening_decision,
    get_screening_statistics
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
"""
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


# Load more functionality removed - now using pagination


@app.route("/privacy")
def privacy():
    return render_template("pandt.html")


@app.route("/support")
def support():
    return render_template("support.html")


@app.route("/search", methods=["GET", "POST"])
def search_page():
    results, source, query, status_filter = [], "", "", None
    page = request.args.get('page', 1, type=int)
    per_page = 5  # Fixed at 5 results per page

    if request.method == "POST":
        source = request.form.get("source")
        query = request.form.get("query", "").strip()
        status_filter = request.form.get("status", "").strip() if source == "clinicaltrials" else None
        page = 1  # Reset to first page for new searches

        if not source or not query:
            flash("⚠ Please select a source and enter a query.", "error")
            return render_template("index.html", results=[], source=source, query=query, pagination={})

        filters = {"status": status_filter} if status_filter else {}
        all_results = search_with_cache(source, query, filters)

        # Sort results alphabetically by title
        all_results.sort(key=lambda x: x.get("title", "").lower())

        # Calculate pagination
        total_results = len(all_results)
        total_pages = (total_results + per_page - 1) // per_page
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page
        results = all_results[start_idx:end_idx]

        # Process results
        for article in results:
            # Standardize PDF field: PMC or OA PDF
            pdf_url = article.get("pmc_pdf_url") or article.get("oa_pdf_path")
            article["oa_pdf_path"] = pdf_url  # overwrite for consistency

            # Format citation
            article["citation"] = format_citation(article)

            # Debug print
            print(f"PDF URL ({source}):", article.get("oa_pdf_path"))

        # Create enhanced pagination info with circular navigation
        pagination = {
            'current_page': page,
            'total_pages': total_pages,
            'total_results': total_results,
            'per_page': per_page,
            'has_prev': page > 1,
            'has_next': page < total_pages,
            'prev_page': page - 1 if page > 1 else total_pages,  # Circular: go to last page if at first
            'next_page': page + 1 if page < total_pages else 1,  # Circular: go to first page if at last
            'start_result': start_idx + 1 if total_results > 0 else 0,
            'end_result': min(end_idx, total_results),
            'page_range': list(range(1, min(total_pages + 1, 101))),  # Show pages 1-100 max
            'show_start': page > 1,
            'is_circular': True  # Enable circular navigation
        }

    elif request.method == "GET" and request.args.get('source') and request.args.get('query'):
        # Handle pagination for existing search results
        source = request.args.get('source')
        query = request.args.get('query')
        status_filter = request.args.get('status')

        filters = {"status": status_filter} if status_filter else {}
        all_results = search_with_cache(source, query, filters)

        # Sort results alphabetically by title
        all_results.sort(key=lambda x: x.get("title", "").lower())

        # Calculate pagination
        total_results = len(all_results)
        total_pages = (total_results + per_page - 1) // per_page
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page
        results = all_results[start_idx:end_idx]

        # Process results
        for article in results:
            # Standardize PDF field: PMC or OA PDF
            pdf_url = article.get("pmc_pdf_url") or article.get("oa_pdf_path")
            article["oa_pdf_path"] = pdf_url  # overwrite for consistency

            # Format citation
            article["citation"] = format_citation(article)

        # Create enhanced pagination info with circular navigation
        pagination = {
            'current_page': page,
            'total_pages': total_pages,
            'total_results': total_results,
            'per_page': per_page,
            'has_prev': page > 1,
            'has_next': page < total_pages,
            'prev_page': page - 1 if page > 1 else total_pages,  # Circular: go to last page if at first
            'next_page': page + 1 if page < total_pages else 1,  # Circular: go to first page if at last
            'start_result': start_idx + 1 if total_results > 0 else 0,
            'end_result': min(end_idx, total_results),
            'page_range': list(range(1, min(total_pages + 1, 101))),  # Show pages 1-100 max
            'show_start': page > 1,
            'is_circular': True  # Enable circular navigation
        }
    else:
        pagination = {}

    return render_template(
        "index.html",
        results=results,
        source=source,
        query=query,
<<<<<<< HEAD
        status=status_filter,
        pagination=pagination
    )


@app.route("/search-history")
def search_history():
    """Enhanced search history with local and global search capabilities."""
    page = request.args.get('page', 1, type=int)
    per_page = 5
    search_type = request.args.get('type', 'local')  # 'local' or 'global'
    search_query = request.args.get('q', '').strip()
    source_filter = request.args.get('source', '')

    # Get filtered history data using enhanced database function
    if search_query or source_filter:
        filtered_history = get_search_history_by_type(
            search_type=search_type,
            search_query=search_query,
            source_filter=source_filter
        )
    else:
        # Get all history with enhanced data structure
        filtered_history = get_archive()

    # Sort alphabetically by query
    filtered_history.sort(key=lambda x: x.get('query', '').lower())

    # Add search type information and statistics
    for item in filtered_history:
        item['display_type'] = search_type
        item['is_local_match'] = search_query and search_query.lower() in item.get('query', '').lower()
        item['is_global_match'] = search_query and (
            search_query.lower() in item.get('query', '').lower() or
            search_query.lower() in item.get('source', '').lower() or
            search_query.lower() in str(item.get('timestamp', '')).lower()
        )

    # Calculate pagination
    total_results = len(filtered_history)
    total_pages = (total_results + per_page - 1) // per_page if total_results > 0 else 1
    start_idx = (page - 1) * per_page
    end_idx = start_idx + per_page
    history_page = filtered_history[start_idx:end_idx]

    # Create enhanced pagination info
    pagination = {
        'current_page': page,
        'total_pages': total_pages,
        'total_results': total_results,
        'per_page': per_page,
        'has_prev': page > 1,
        'has_next': page < total_pages,
        'prev_page': page - 1 if page > 1 else total_pages,  # Circular navigation
        'next_page': page + 1 if page < total_pages else 1,  # Circular navigation
        'start_result': start_idx + 1 if total_results > 0 else 0,
        'end_result': min(end_idx, total_results),
        'page_range': list(range(1, min(total_pages + 1, 101))),  # Show pages 1-100 max
        'show_start': page > 1,
        'is_circular': True
    }

    # Get available sources for filter
    available_sources = list(set(h.get('source', '') for h in all_history if h.get('source')))
    available_sources.sort()

    return render_template("search_history.html",
                         history=history_page,
                         pagination=pagination,
                         search_query=search_query,
                         search_type=search_type,
                         source_filter=source_filter,
                         available_sources=available_sources)
=======
        page=1,             # current page number (default to 1)
        per_page=20,        # results per page
        total=len(results), # total results count (use API total if available)
        status=status_filter
    )



@app.route("/history")
def history():
    history_data = get_archive()
    return render_template("history.html", history=history_data)
>>>>>>> d3bf0ae (Your update)


@app.route("/donations")
def donations():
    return render_template("donations.html")


@app.route("/blogs")
def blogs():
    return render_template("blogs.html")

<<<<<<< HEAD
@app.route('/home')
def home():
    return render_template('index.html')  # Or your main search page template
=======
@app.route('/home') #This section updates for database error
def home():
    return render_template(
        'index.html',# Or your main search page template
        page=1,
        results=[],
        total=0,
        query=""
    )
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
   # return {
    #    'results': id_list,
     #   'page': page,
      #  'per_page': per_page,
       # 'total': total_count
   # }

    # Call the actual search function (you need to ensure this function exists and accepts these params)
    result = search_pubmed(query, page, per_page)

    print("Returning result:", result)

    return jsonify(result)
>>>>>>> d3bf0ae (Your update)


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


# Pagination and search routes removed as requested


# -------------------------------------------------
# Screening Project Management Routes
# -------------------------------------------------
@app.route("/screening")
def screening_dashboard():
    """Main screening dashboard showing all projects."""
    from database.database import get_screening_projects
    projects = get_screening_projects()
    return render_template("screening/dashboard.html", projects=projects)


@app.route("/screening/create", methods=["GET", "POST"])
def create_screening_project():
    """Create a new screening project."""
    if request.method == "POST":
        from database.database import create_screening_project, add_reviewer_to_project

        name = request.form.get("name")
        description = request.form.get("description", "")
        inclusion_criteria = request.form.get("inclusion_criteria", "")
        exclusion_criteria = request.form.get("exclusion_criteria", "")
        created_by = request.form.get("created_by", "Anonymous")
        require_dual_screening = bool(request.form.get("require_dual_screening"))
        conflict_resolution_method = request.form.get("conflict_resolution_method", "discussion")

        # Create project
        project_id = create_screening_project(
            name, description, inclusion_criteria, exclusion_criteria,
            created_by, require_dual_screening, conflict_resolution_method
        )

        if project_id:
            # Add creator as lead reviewer
            add_reviewer_to_project(project_id, created_by, "", "lead_reviewer")

            # Check if we should import search results
            import_search_results = request.form.get("import_search_results")
            if import_search_results == "true":
                search_source = request.form.get("search_source")
                search_query = request.form.get("search_query")

                # Import current search results into the project
                if search_source and search_query:
                    try:
                        # Get the most recent search results for this query
                        recent_search = None
                        search_history = get_archive()
                        for search in search_history:
                            if search.source == search_source and search.query == search_query:
                                recent_search = search
                                break

                        if recent_search:
                            articles = get_results_by_search_id(recent_search.id)
                            # Convert to screening format
                            screening_articles = []
                            for article in articles:
                                screening_articles.append({
                                    'title': article[0],
                                    'authors': article[1],
                                    'abstract': article[5],
                                    'doi': article[2],
                                    'pmid': article[3],
                                    'url': article[4]
                                })

                            count = add_articles_to_screening_project(project_id, screening_articles)
                            flash(f"✅ Screening project '{name}' created with {count} articles imported!", "success")
                        else:
                            flash(f"✅ Screening project '{name}' created successfully! No recent search results found to import.", "info")
                    except Exception as e:
                        print(f"Error importing search results: {e}")
                        flash(f"✅ Screening project '{name}' created successfully! Error importing search results.", "warning")
                else:
                    flash(f"✅ Screening project '{name}' created successfully!", "success")
            else:
                flash(f"✅ Screening project '{name}' created successfully!", "success")

            return redirect(url_for("screening_project", project_id=project_id))
        else:
            flash("❌ Failed to create screening project.", "error")

    return render_template("screening/create_project.html")


@app.route("/screening/project/<int:project_id>")
def screening_project(project_id):
    """View a specific screening project."""
    from database.database import (get_screening_project, get_project_reviewers,
                                   get_screening_statistics, get_screening_articles)

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    reviewers = get_project_reviewers(project_id)
    statistics = get_screening_statistics(project_id)
    articles = get_screening_articles(project_id)

    return render_template("screening/project_detail.html",
                         project=project, reviewers=reviewers,
                         statistics=statistics, articles=articles)


@app.route("/screening/project/<int:project_id>/add_reviewer", methods=["POST"])
def add_reviewer(project_id):
    """Add a reviewer to a screening project."""
    from database.database import add_reviewer_to_project

    reviewer_name = request.form.get("reviewer_name")
    reviewer_email = request.form.get("reviewer_email", "")
    role = request.form.get("role", "reviewer")

    if add_reviewer_to_project(project_id, reviewer_name, reviewer_email, role):
        flash(f"✅ Reviewer '{reviewer_name}' added successfully!", "success")
    else:
        flash("❌ Failed to add reviewer.", "error")

    return redirect(url_for("screening_project", project_id=project_id))


@app.route("/screening/project/<int:project_id>/add_articles", methods=["GET", "POST"])
def add_articles_to_screening(project_id):
    """Add articles to a screening project from search results."""
    from database.database import get_screening_project, add_articles_to_screening_project

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    if request.method == "POST":
        # Handle adding articles from search results or manual input
        source = request.form.get("source")

        if source == "search_results":
            # Get articles from a previous search
            search_id = request.form.get("search_id")
            if search_id:
                articles = get_results_by_search_id(search_id)
                # Convert to screening format
                screening_articles = []
                for article in articles:
                    screening_articles.append({
                        'title': article[0],
                        'authors': article[1],
                        'abstract': article[5],
                        'doi': article[2],
                        'pmid': article[3],
                        'url': article[4]
                    })

                count = add_articles_to_screening_project(project_id, screening_articles)
                flash(f"✅ Added {count} articles to screening project!", "success")
                return redirect(url_for("screening_project", project_id=project_id))

        elif source == "manual":
            # Handle manual article input
            title = request.form.get("title")
            authors = request.form.get("authors", "")
            abstract = request.form.get("abstract", "")
            doi = request.form.get("doi", "")
            pmid = request.form.get("pmid", "")
            url = request.form.get("url", "")

            if title:
                articles = [{
                    'title': title,
                    'authors': authors,
                    'abstract': abstract,
                    'doi': doi,
                    'pmid': pmid,
                    'url': url
                }]

                count = add_articles_to_screening_project(project_id, articles)
                flash(f"✅ Added {count} article to screening project!", "success")
                return redirect(url_for("screening_project", project_id=project_id))

    # Get available search history for adding articles
    search_history = get_archive()
    return render_template("screening/add_articles.html",
                         project=project, search_history=search_history)


@app.route("/screening/project/<int:project_id>/screen/<int:article_id>", methods=["GET", "POST"])
def screen_article(project_id, article_id):
    """Screen an individual article."""
    from database.database import (get_screening_project, get_screening_articles,
                                   save_screening_decision, get_project_reviewers)

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    # Get the specific article
    articles = get_screening_articles(project_id)
    article = None
    for a in articles:
        if a.id == article_id:
            article = a
            break

    if not article:
        flash("❌ Article not found.", "error")
        return redirect(url_for("screening_project", project_id=project_id))

    reviewers = get_project_reviewers(project_id)

    if request.method == "POST":
        reviewer_id = request.form.get("reviewer_id")
        decision = request.form.get("decision")
        reason = request.form.get("reason", "")
        notes = request.form.get("notes", "")
        screening_stage = request.form.get("screening_stage", "title_abstract")

        if save_screening_decision(article_id, reviewer_id, decision, reason, notes, screening_stage):
            flash(f"✅ Decision '{decision}' saved successfully!", "success")

            # Redirect to next article or back to project
            next_action = request.form.get("next_action", "project")
            if next_action == "next_article":
                # Find next unscreened article
                unscreened = get_screening_articles(project_id, status="pending", reviewer_id=reviewer_id)
                if unscreened and len(unscreened) > 0:
                    return redirect(url_for("screen_article", project_id=project_id, article_id=unscreened[0].id))

            return redirect(url_for("screening_project", project_id=project_id))
        else:
            flash("❌ Failed to save decision.", "error")

    return render_template("screening/screen_article.html",
                         project=project, article=article, reviewers=reviewers)


@app.route("/screening/project/<int:project_id>/batch_screen", methods=["GET", "POST"])
def batch_screen(project_id):
    """Batch screening interface for efficient screening."""
    from database.database import (get_screening_project, get_screening_articles,
                                   save_screening_decision, get_project_reviewers)

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    reviewers = get_project_reviewers(project_id)

    if request.method == "POST":
        reviewer_id = request.form.get("reviewer_id")
        decisions = request.form.getlist("decisions")
        article_ids = request.form.getlist("article_ids")

        saved_count = 0
        for i, article_id in enumerate(article_ids):
            if i < len(decisions) and decisions[i]:
                decision_data = decisions[i].split("|")
                if len(decision_data) >= 2:
                    decision = decision_data[0]
                    reason = decision_data[1] if len(decision_data) > 1 else ""

                    if save_screening_decision(article_id, reviewer_id, decision, reason):
                        saved_count += 1

        flash(f"✅ Saved {saved_count} screening decisions!", "success")
        return redirect(url_for("screening_project", project_id=project_id))

    # Get pending articles for batch screening
    articles = get_screening_articles(project_id, status="pending")

    return render_template("screening/batch_screen.html",
                         project=project, articles=articles, reviewers=reviewers)


@app.route("/screening/project/<int:project_id>/conflicts")
def view_conflicts(project_id):
    """View and resolve screening conflicts."""
    from database.database import get_screening_project

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    # Get articles with conflicts
    conflict_articles = get_screening_articles(project_id, status="conflict")

    return render_template("screening/conflicts.html",
                         project=project, conflict_articles=conflict_articles)


@app.route("/screening/project/<int:project_id>/analytics")
def screening_analytics(project_id):
    """View detailed analytics for a screening project."""
    from database.database import (get_screening_project, get_screening_statistics,
                                   get_project_reviewers, get_screening_articles)

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    statistics = get_screening_statistics(project_id)
    reviewers = get_project_reviewers(project_id)

    # Calculate additional analytics
    analytics_data = calculate_screening_analytics(project_id)

    return render_template("screening/analytics.html",
                         project=project, statistics=statistics,
                         reviewers=reviewers, analytics=analytics_data)


@app.route("/screening/project/<int:project_id>/export")
def export_screening_results(project_id):
    """Export screening results in various formats."""
    from database.database import get_screening_project, get_screening_articles
    import csv
    import io
    from flask import make_response

    project = get_screening_project(project_id)
    if not project:
        flash("❌ Screening project not found.", "error")
        return redirect(url_for("screening_dashboard"))

    export_format = request.args.get('format', 'csv')
    status_filter = request.args.get('status', 'all')

    # Get articles based on filter
    if status_filter == 'all':
        articles = get_screening_articles(project_id)
    else:
        articles = get_screening_articles(project_id, status=status_filter)

    if export_format == 'csv':
        output = io.StringIO()
        writer = csv.writer(output)

        # Write header
        writer.writerow(['Title', 'Authors', 'DOI', 'PMID', 'Abstract', 'Status', 'URL'])

        # Write data
        for article in articles:
            writer.writerow([
                article.title,
                article.authors,
                article.doi or '',
                article.pmid or '',
                article.abstract or '',
                article.screening_status,
                article.url or ''
            ])

        output.seek(0)
        response = make_response(output.getvalue())
        response.headers['Content-Type'] = 'text/csv'
        response.headers['Content-Disposition'] = f'attachment; filename={project.name}_screening_results.csv'
        return response

    # For other formats, redirect back with message
    flash("Export format not yet supported. CSV export is available.", "info")
    return redirect(url_for("screening_project", project_id=project_id))


def calculate_screening_analytics(project_id):
    """Calculate detailed analytics for screening project."""
    from database.database import get_connection

    analytics = {}

    try:
        with get_connection() as conn:
            cursor = conn.cursor()

            # Inter-rater agreement calculation
            cursor.execute("""
                SELECT sa.id, sd.reviewer_id, sd.decision
                FROM screening_articles sa
                JOIN screening_decisions sd ON sa.id = sd.screening_article_id
                WHERE sa.project_id = ?
                ORDER BY sa.id, sd.reviewer_id
            """, (project_id,))

            decisions_by_article = {}
            for row in cursor.fetchall():
                article_id, reviewer_id, decision = row
                if article_id not in decisions_by_article:
                    decisions_by_article[article_id] = {}
                decisions_by_article[article_id][reviewer_id] = decision

            # Calculate agreement statistics
            total_agreements = 0
            total_comparisons = 0

            for article_id, decisions in decisions_by_article.items():
                if len(decisions) >= 2:
                    decision_values = list(decisions.values())
                    for i in range(len(decision_values)):
                        for j in range(i + 1, len(decision_values)):
                            total_comparisons += 1
                            if decision_values[i] == decision_values[j]:
                                total_agreements += 1

            agreement_rate = (total_agreements / total_comparisons * 100) if total_comparisons > 0 else 0
            analytics['agreement_rate'] = round(agreement_rate, 1)
            analytics['total_comparisons'] = total_comparisons

            # Screening velocity (decisions per day)
            cursor.execute("""
                SELECT DATE(sd.decision_date) as decision_day, COUNT(*) as decisions_count
                FROM screening_decisions sd
                JOIN screening_articles sa ON sd.screening_article_id = sa.id
                WHERE sa.project_id = ?
                GROUP BY DATE(sd.decision_date)
                ORDER BY decision_day DESC
                LIMIT 7
            """, (project_id,))

            daily_decisions = cursor.fetchall()
            analytics['daily_decisions'] = daily_decisions

            # Reviewer performance
            cursor.execute("""
                SELECT sr.reviewer_name,
                       COUNT(sd.id) as total_decisions,
                       COUNT(CASE WHEN sd.decision = 'include' THEN 1 END) as include_count,
                       COUNT(CASE WHEN sd.decision = 'exclude' THEN 1 END) as exclude_count,
                       COUNT(CASE WHEN sd.decision = 'maybe' THEN 1 END) as maybe_count
                FROM screening_reviewers sr
                LEFT JOIN screening_decisions sd ON sr.id = sd.reviewer_id
                WHERE sr.project_id = ?
                GROUP BY sr.id, sr.reviewer_name
            """, (project_id,))

            reviewer_performance = cursor.fetchall()
            analytics['reviewer_performance'] = reviewer_performance

    except Exception as e:
        print(f"[Analytics] ⚠ Error calculating analytics: {e}")
        analytics = {
            'agreement_rate': 0,
            'total_comparisons': 0,
            'daily_decisions': [],
            'reviewer_performance': []
        }

    return analytics


# -------------------------------------------------
# Main Entry
# -------------------------------------------------
if __name__ == "__main__":
    init_db()
    app.run(host="0.0.0.0", port=5000, debug=True)



