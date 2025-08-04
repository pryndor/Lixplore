# database.py
import sqlite3
import os
import json
from datetime import datetime

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

        # ✅ Enable WAL mode for better performance
        cursor.execute("PRAGMA journal_mode=WAL;")

        # ✅ Create Searches Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS searches (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source TEXT NOT NULL,            
            query TEXT NOT NULL,
            filters TEXT NOT NULL,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
        )
        """)

        # ✅ Create Articles Table
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
    """Return SQLite connection with WAL mode enabled."""
    conn = sqlite3.connect(DB_PATH, timeout=10)
    conn.execute("PRAGMA journal_mode=WAL;")
    return conn


# -------------------------------------------------
# Save Search & Articles
# -------------------------------------------------
def save_search(source, query, filters, results):
    """
    Save search query and related articles to the database.
    Filters are stored as JSON string with sorted keys for consistency.
    """
    try:
        filters_str = json.dumps(filters or {}, sort_keys=True)

        with get_connection() as conn:
            cursor = conn.cursor()

            # ✅ Insert search entry
            cursor.execute("""
                INSERT INTO searches (source, query, filters) 
                VALUES (?, ?, ?)
            """, (source, query, filters_str))
            search_id = cursor.lastrowid
            print(f"[DB] ✅ Saved search ID={search_id}: {source} | {query}")

            # ✅ Insert all articles
            for article in results:
                if not isinstance(article, dict):
                    print(f"[DB] ⚠ Skipping invalid article (not a dict): {article}")
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
            print(f"[DB] ✅ Saved {len(results)} articles for search ID={search_id}")

    except sqlite3.Error as e:
        print(f"[DB] ⚠ SQLite error while saving search: {e}")
    except Exception as e:
        print(f"[DB] ⚠ Unexpected error while saving search: {e}")


# -------------------------------------------------
# Retrieve Search History
# -------------------------------------------------
def get_archive():
    """Fetch all saved search history in descending order of timestamp."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT id, source, query, filters, timestamp 
                FROM searches 
                ORDER BY timestamp DESC
            """)
            history = cursor.fetchall()
            print(f"[DB] 📜 Retrieved {len(history)} search history entries")
            return history
    except sqlite3.Error as e:
        print(f"[DB] ⚠ Failed to retrieve search history: {e}")
        return []


# -------------------------------------------------
# Retrieve Results by Search ID
# -------------------------------------------------
def get_results_by_search_id(search_id):
    """Fetch all articles related to a saved search ID."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                FROM articles 
                WHERE search_id=?
            """, (search_id,))
            rows = cursor.fetchall()
            print(f"[DB] 📜 Retrieved {len(rows)} articles for search ID={search_id}")
            return rows
    except sqlite3.Error as e:
        print(f"[DB] ⚠ Failed to retrieve articles for search {search_id}: {e}")
        return []

