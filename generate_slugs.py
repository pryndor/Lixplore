import sqlite3
import re

def slugify(title):
    slug = re.sub(r'[\W_]+', '-', title.lower()).strip('-')
    return slug

DB_PATH = "database/lixplore.db"

with sqlite3.connect(DB_PATH) as conn:
    cursor = conn.cursor()
    cursor.execute("SELECT id, title FROM blog_posts WHERE slug IS NULL OR slug = ''")
    posts = cursor.fetchall()

    for post_id, title in posts:
        slug = slugify(title)
        print(f"Generating slug '{slug}' for post ID {post_id}")
        cursor.execute("UPDATE blog_posts SET slug = ? WHERE id = ?", (slug, post_id))

    conn.commit()
    print("✅ All missing slugs updated.")

