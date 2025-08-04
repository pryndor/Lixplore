# diagnostic_check.py
import sqlite3
import os

# Path to your DB
DB_PATH = os.path.join(os.path.dirname(__file__), "literature_archive.db")

def check_article_counts():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT s.id, s.source, s.query, COUNT(a.id) as article_count
        FROM searches s
        LEFT JOIN articles a ON s.id = a.search_id
        GROUP BY s.id
        ORDER BY s.timestamp DESC
    """)

    rows = cursor.fetchall()

    print("\n=== Search History Article Count Diagnostic ===")
    for row in rows:
        search_id, source, query, article_count = row
        status = "✅ OK" if article_count > 0 else "❌ NO ARTICLES"
        print(f"[{status}] ID={search_id} | {source} | {query} | Articles={article_count}")

    conn.close()

if __name__ == "__main__":
    check_article_counts()

