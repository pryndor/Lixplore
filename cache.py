# cache.py
import json
import sqlite3
from datetime import datetime, timedelta

DB_PATH = "/home/bala/Lixplore/database/lixplore.db"

def get_connection():
    """Create and return a new database connection."""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def get_cached_results(source, query, filters):
    """
    Check cache for existing results for given source, query, and filters.
    Returns a tuple: (results as list/dict, search_cache id) or (None, None) if no valid cache.
    """
    filters_json = json.dumps(filters, sort_keys=True) if filters else None

    conn = get_connection()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT id, results, created_at FROM search_cache
        WHERE source = ? AND query = ? AND (filters IS ? OR filters = ?)
        ORDER BY created_at DESC LIMIT 1
    """, (source, query, filters_json, filters_json))
    row = cursor.fetchone()
    conn.close()

    if row:
        created_at = datetime.strptime(row["created_at"], "%Y-%m-%d %H:%M:%S")
        # Cache expiry: 7 days (change as needed)
        if datetime.now() - created_at > timedelta(days=7):
            return None, None  # Cache expired

        return json.loads(row["results"]), row["id"]

    return None, None  # No cache found

def save_results_to_cache(source, query, filters, results):
    """
    Save search results JSON to cache and save detailed articles.
    Returns the inserted search_cache record ID.
    """
    filters_json = json.dumps(filters, sort_keys=True) if filters else None
    results_json = json.dumps(results)

    conn = get_connection()
    cursor = conn.cursor()
    cursor.execute("""
        INSERT INTO search_cache (source, query, filters, results, created_at)
        VALUES (?, ?, ?, ?, CURRENT_TIMESTAMP)
    """, (source, query, filters_json, results_json))
    search_id = cursor.lastrowid

    # Save detailed articles related to this search
    for article in results:
        cursor.execute("""
            INSERT INTO articles (
                search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            search_id,
            article.get("title"),
            ", ".join(article.get("authors", [])) if article.get("authors") else None,
            article.get("doi"),
            article.get("pmid"),
            article.get("url"),
            article.get("abstract"),
            ", ".join(article.get("affiliations", [])) if article.get("affiliations") else None,
            article.get("oa_pdf_path"),
        ))

    conn.commit()
    conn.close()

    return search_id

