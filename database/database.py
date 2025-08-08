

import sqlite3
import os
import json
from datetime import datetime


# 📂 Database file path
BASE_DIR = os.path.dirname(__file__)  # This points to /home/bala/Lixplore/database
DB_PATH = os.path.join(BASE_DIR, "lixplore.db")  # Correctly points to /home/bala/Lixplore/database/lixplore.db


# -------------------------------------------------
# Initialize Database
# -------------------------------------------------
def init_db():
    """Initialize database with WAL mode and ensure schema is correct."""
    with sqlite3.connect(DB_PATH, timeout=10) as conn:
        cursor = conn.cursor()
        cursor.execute("PRAGMA journal_mode=WAL;")

        # Searches Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS searches (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source TEXT NOT NULL,            
            query TEXT NOT NULL,
            filters TEXT NOT NULL,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
        )
        """)

        # Articles Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS articles (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            search_id INTEGER NOT NULL,
            title TEXT,
            authors TEXT,
            doi TEXT,
            pmid TEXT,
            url TEXT,
            abstract TEXT,
            affiliations TEXT,
            oa_pdf_path TEXT,
            FOREIGN KEY (search_id) REFERENCES searches(id) ON DELETE CASCADE
        )
        """)

        # Blogs Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS blogs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            slug TEXT UNIQUE NOT NULL,
            title TEXT NOT NULL,
            author TEXT,
            content TEXT NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )
        """)

        # Blog Tags Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS blog_tags (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE NOT NULL
        )
        """)

        # Blog-Tag Mapping Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS blog_tag_map (
            blog_id INTEGER,
            tag_id INTEGER,
            PRIMARY KEY (blog_id, tag_id),
            FOREIGN KEY (blog_id) REFERENCES blogs(id) ON DELETE CASCADE,
            FOREIGN KEY (tag_id) REFERENCES blog_tags(id) ON DELETE CASCADE
        )
        """)

        conn.commit()
        print("[DB] ✅ Database initialized successfully")


# -------------------------------------------------
# Connection Helper
# -------------------------------------------------
def get_connection():
    conn = sqlite3.connect(DB_PATH, timeout=10)
    conn.row_factory = sqlite3.Row
    return conn


# -------------------------------------------------
# Blog CRUD
# -------------------------------------------------
def save_blog(slug, title, author, content, tags=None):
    """Save a blog post with optional tags."""
    tags = tags or []

    with get_connection() as conn:
        cursor = conn.cursor()

        # Insert blog post
        cursor.execute("""
            INSERT INTO blogs (slug, title, author, content, created_at)
            VALUES (?, ?, ?, ?, ?)
        """, (slug, title, author, content, datetime.utcnow()))

        blog_id = cursor.lastrowid

        # Insert tags and mapping
        for tag in tags:
            # Ensure tag exists or insert it
            cursor.execute("INSERT OR IGNORE INTO blog_tags (name) VALUES (?)", (tag,))
            cursor.execute("SELECT id FROM blog_tags WHERE name = ?", (tag,))
            tag_id = cursor.fetchone()[0]

            # Map blog to tag
            cursor.execute("""
                INSERT OR IGNORE INTO blog_tag_map (blog_id, tag_id) VALUES (?, ?)
            """, (blog_id, tag_id))

        conn.commit()
        print(f"[DB] ✅ Blog saved: {title} | Tags: {tags}")

def get_posts_by_tag(tag):
    conn = get_connection()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT * FROM blog_posts
        WHERE tags LIKE ?
        ORDER BY created_at DESC
    """, (f"%{tag}%",))
    posts = cursor.fetchall()
    conn.close()
    return posts



def get_blog_by_slug(slug):
    """Fetch blog post by slug including tags."""
    with get_connection() as conn:
        cursor = conn.cursor()

        # Fetch blog
        cursor.execute("""
            SELECT id, title, author, content, created_at FROM blogs WHERE slug = ?
        """, (slug,))
        row = cursor.fetchone()
        if not row:
            return None

        blog_id, title, author, content, created_at = row

        # Fetch tags
        cursor.execute("""
            SELECT name FROM blog_tags 
            JOIN blog_tag_map ON blog_tags.id = blog_tag_map.tag_id
            WHERE blog_tag_map.blog_id = ?
        """, (blog_id,))
        tags = [tag_row[0] for tag_row in cursor.fetchall()]

        return {
            "title": title,
            "author": author,
            "content": content,
            "created_at": created_at,
            "tags": tags
        }


def get_all_blogs():
    """Fetch all blog posts with tags."""
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT id, slug, title, author, created_at FROM blogs ORDER BY created_at DESC
        """)
        blogs = []
        for row in cursor.fetchall():
            blog_id, slug, title, author, created_at = row
            cursor.execute("""
                SELECT name FROM blog_tags 
                JOIN blog_tag_map ON blog_tags.id = blog_tag_map.tag_id
                WHERE blog_tag_map.blog_id = ?
            """, (blog_id,))
            tags = [tag_row[0] for tag_row in cursor.fetchall()]
            blogs.append({
                "slug": slug,
                "title": title,
                "author": author,
                "created_at": created_at,
                "tags": tags
            })
        return blogs


# -------------------------------------------------
# Existing search-related functions (unchanged)
# -------------------------------------------------

def save_search(source, query, filters, results):
    try:
        filters_str = json.dumps(filters or {}, sort_keys=True)
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO searches (source, query, filters) 
                VALUES (?, ?, ?)
            """, (source, query, filters_str))
            search_id = cursor.lastrowid

            for article in (results or []):
                if not isinstance(article, dict):
                    continue
                cursor.execute("""
                    INSERT INTO articles 
                    (search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    search_id,
                    article.get("title") or "Untitled",
                    ", ".join(article.get("authors", [])) if isinstance(article.get("authors"), list) else (article.get("authors") or ""),
                    article.get("doi"),
                    article.get("pmid"),
                    article.get("url"),
                    article.get("abstract"),
                    ", ".join(article.get("affiliations", [])) if isinstance(article.get("affiliations"), list) else (article.get("affiliations") or ""),
                    article.get("oa_pdf_path")
                ))

            conn.commit()
            print(f"[DB] ✅ Articles saved for search ID={search_id}")

    except Exception as e:
        print(f"[DB] ⚠ Error saving search: {e}")


def get_archive():
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT id, source, query, filters, timestamp 
                FROM searches 
                ORDER BY timestamp DESC
            """)
            return cursor.fetchall()
    except Exception as e:
        print(f"[DB] ⚠ Failed to retrieve history: {e}")
        return []


def get_results_by_search_id(search_id):
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                FROM articles 
                WHERE search_id=?
            """, (search_id,))
            return cursor.fetchall()
    except Exception as e:
        print(f"[DB] ⚠ Failed to retrieve results: {e}")
        return []

