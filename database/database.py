# Sqllite
import sqlite3
import os

# 📂 Database file path
DB_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "literature_archive.db")


def init_db():
    """Initialize database with WAL mode for better concurrency."""
    with sqlite3.connect(DB_PATH, timeout=10) as conn:
        cursor = conn.cursor()

        # ✅ Enable WAL mode (better concurrency, fewer lock issues)
        cursor.execute("PRAGMA journal_mode=WAL;")

        # ✅ Create Searches Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS searches (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source TEXT,            
            query TEXT,
            filters TEXT,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
        )
        """)

        # ✅ Create Articles Table
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS articles (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            search_id INTEGER,
            title TEXT,
            authors TEXT,
            doi TEXT,
            pmid TEXT,
            url TEXT,
            abstract TEXT,
            affiliations TEXT,
            oa_pdf_path TEXT,       
            FOREIGN KEY (search_id) REFERENCES searches(id)
        )
        """)

        conn.commit()


def get_connection():
    """Return SQLite connection with WAL mode and timeout."""
    conn = sqlite3.connect(DB_PATH, timeout=10)
    conn.execute("PRAGMA journal_mode=WAL;")
    return conn


def save_search(source, query, filters, results):
    """Save search query and related articles to the database."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()

            # ✅ Insert into searches table
            cursor.execute(
                "INSERT INTO searches (source, query, filters) VALUES (?, ?, ?)",
                (source, query, str(filters))
            )
            search_id = cursor.lastrowid

            # ✅ Insert related articles
            for article in results:
                # Defensive check
                if not isinstance(article, dict):
                    print(f"⚠ Skipping invalid article (not dict): {article}")
                    continue  

                cursor.execute("""
                INSERT INTO articles 
                (search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    search_id,
                    article.get("title"),
                    ", ".join(article.get("authors", [])) if isinstance(article.get("authors"), list) else article.get("authors"),
                    article.get("doi"),
                    article.get("pmid"),
                    article.get("url"),
                    article.get("abstract"),
                    ", ".join(article.get("affiliations", [])) if isinstance(article.get("affiliations"), list) else article.get("affiliations"),
                    article.get("oa_pdf_path")
                ))

            conn.commit()

    except sqlite3.OperationalError as e:
        print(f"⚠ Database save error: {e}")
    except Exception as e:
        print(f"⚠ Unexpected error while saving search: {e}")


def get_archive():
    """Fetch all saved search history."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT id, source, query, filters, timestamp 
                FROM searches 
                ORDER BY timestamp DESC
            """)
            return cursor.fetchall()
    except sqlite3.Error as e:
        print(f"⚠ Failed to retrieve search history: {e}")
        return []


