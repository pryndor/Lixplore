# blog.py

import sqlite3
from flask import Blueprint, render_template, abort
from markdown import markdown

blog_bp = Blueprint('blog', __name__)

DB_PATH = "database/literature_archive.db"

def get_blog_by_slug(slug):
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute("SELECT title, author, content, created_at FROM blogs WHERE slug = ?", (slug,))
    row = cursor.fetchone()
    conn.close()
    if row:
        return {
            "title": row[0],
            "author": row[1],
            "content": markdown(row[2]),  # Convert Markdown to HTML
            "created_at": row[3]
        }
    return None

def get_all_blogs():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute("SELECT slug, title, author, created_at FROM blogs ORDER BY created_at DESC")
    rows = cursor.fetchall()
    conn.close()
    return rows

@blog_bp.route("/blog")
def blog_list():
    blogs = get_all_blogs()
    return render_template("blog_list.html", blogs=blogs)

@blog_bp.route("/blog/<slug>")
def show_blog(slug):
    blog = get_blog_by_slug(slug)
    if not blog:
        abort(404)
    return render_template("blog.html", blog=blog)

