#!/usr/bin/env python3
import sqlite3
import os
import markdown
import re
import yaml
from datetime import datetime

DB_PATH = "database/lixplore.db"
BLOG_DIR = "blog_posts/"

FRONTMATTER_RE = re.compile(r'^---\s*\n(.*?)\n---\s*\n(.*)$', re.DOTALL)

def slugify(title):
    return re.sub(r'[\W_]+', '-', title.lower()).strip('-')

def parse_markdown_file(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    m = FRONTMATTER_RE.match(content)
    if m:
        fm_text, body = m.group(1), m.group(2)
        try:
            frontmatter = yaml.safe_load(fm_text) or {}
        except Exception:
            frontmatter = {}
    else:
        frontmatter = {}
        body = content

    html = markdown.markdown(body)
    return frontmatter, html

def insert_blog_from_md():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    for filename in sorted(os.listdir(BLOG_DIR)):
        if not filename.endswith(".md"):
            continue
        filepath = os.path.join(BLOG_DIR, filename)
        frontmatter, html = parse_markdown_file(filepath)

        title = frontmatter.get("title") or filename.replace(".md", "").replace("-", " ").title()
        tags = frontmatter.get("tags", [])
        if isinstance(tags, str):
            tags = [t.strip() for t in tags.split(",") if t.strip()]
        tags_str = ",".join(tags)

        author = frontmatter.get("author", "Admin")
        slug = frontmatter.get("slug") or slugify(title)

        try:
            cursor.execute("SELECT id FROM blog_posts WHERE slug = ?", (slug,))
            row = cursor.fetchone()
            if row:
                blog_id = row[0]
                cursor.execute("""
                    UPDATE blog_posts
                    SET title = ?, author = ?, content = ?, tags = ?
                    WHERE id = ?
                """, (title, author, html, tags_str, blog_id))
                print(f"♻️ Updated: {filename} → /blog/{slug}")
            else:
                cursor.execute("""
                    INSERT INTO blog_posts (title, author, content, tags, slug)
                    VALUES (?, ?, ?, ?, ?)
                """, (title, author, html, tags_str, slug))
                print(f"✅ Inserted: {filename} → /blog/{slug}")
        except Exception as e:
            print(f"❌ Error with {filename}: {e}")

    conn.commit()
    conn.close()

if __name__ == "__main__":
    insert_blog_from_md()

