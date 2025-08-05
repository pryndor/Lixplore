# database.py

import sqlite3
import os
import json

# 📂 Database file path
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DB_PATH = os.path.join(os.path.dirname(__file__), "literature_archive.db")


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

        conn.commit()
        print("[DB] ✅ Database initialized successfully")


# -------------------------------------------------
# Connection Helper
# -------------------------------------------------
def get_connection():
    conn = sqlite3.connect(DB_PATH, timeout=10)
    conn.execute("PRAGMA journal_mode=WAL;")
    return conn


# -------------------------------------------------
# Save Search & Articles
# -------------------------------------------------
def save_search(source, query, filters, results):
    """Save search query & results (even if results are empty)."""
    try:
        filters_str = json.dumps(filters or {}, sort_keys=True)

        with get_connection() as conn:
            cursor = conn.cursor()

            # Save search always
            cursor.execute("""
                INSERT INTO searches (source, query, filters) 
                VALUES (?, ?, ?)
            """, (source, query, filters_str))
            search_id = cursor.lastrowid
            print(f"[DB] ✅ Search saved ID={search_id}: {source} | {query}")

            # Save articles (only if present)
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


# -------------------------------------------------
# Retrieve Search History
# -------------------------------------------------
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


# -------------------------------------------------
# Retrieve Results by Search ID
# -------------------------------------------------
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

