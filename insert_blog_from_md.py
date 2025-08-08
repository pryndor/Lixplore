import sqlite3
import os
import markdown
import re
import yaml

DB_PATH = "database/lixplore.db"
BLOG_DIR = "blog_posts/"

def slugify(title):
    return re.sub(r'[\W_]+', '-', title.lower()).strip('-')

def parse_markdown_file(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    if content.startswith('---'):
        _, fm_text, body = content.split('---', 2)
        frontmatter = yaml.safe_load(fm_text)
    else:
        frontmatter = {}
        body = content

    html = markdown.markdown(body)
    return frontmatter, html

def insert_blog_from_md():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    for filename in os.listdir(BLOG_DIR):
        if filename.endswith(".md"):
            filepath = os.path.join(BLOG_DIR, filename)
            frontmatter, html = parse_markdown_file(filepath)

            title = frontmatter.get("title", "Untitled")
            tags = frontmatter.get("tags", [])
            if isinstance(tags, str):
                tags = [tag.strip() for tag in tags.split(',')]
            tags_str = ",".join(tags)

            author = frontmatter.get("author", "Admin")
            slug = frontmatter.get("slug", slugify(title))

            try:
                cursor.execute("SELECT 1 FROM blog_posts WHERE slug = ?", (slug,))
                if cursor.fetchone():
                    print(f"⚠️  Skipping already indexed: {filename} → /blog/{slug}")
                    continue

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

