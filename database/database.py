

import sqlite3
import os
import json
from datetime import datetime

# 📂 Database file path
BASE_DIR = os.path.dirname(__file__)  # e.g., /home/bala/Lixplore/database
DB_PATH = os.path.join(BASE_DIR, "lixplore.db")  # /home/bala/Lixplore/database/lixplore.db


# -------------------------------------------------
# Initialize Database
# -------------------------------------------------
def init_db():
    """Initialize database in WAL mode and ensure schema is correct."""
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

        # Screening Projects Table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_projects (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                description TEXT,
                inclusion_criteria TEXT,
                exclusion_criteria TEXT,
                created_by TEXT NOT NULL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                status TEXT DEFAULT 'active' CHECK (status IN ('active', 'completed', 'paused')),
                require_dual_screening BOOLEAN DEFAULT 0,
                conflict_resolution_method TEXT DEFAULT 'discussion' CHECK (conflict_resolution_method IN ('discussion', 'third_reviewer', 'consensus'))
            )
        """)

        # Screening Reviewers Table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_reviewers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                project_id INTEGER NOT NULL,
                reviewer_name TEXT NOT NULL,
                reviewer_email TEXT,
                role TEXT DEFAULT 'reviewer' CHECK (role IN ('reviewer', 'lead_reviewer', 'admin')),
                added_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (project_id) REFERENCES screening_projects(id) ON DELETE CASCADE,
                UNIQUE(project_id, reviewer_email)
            )
        """)

        # Screening Articles Table (links articles to screening projects)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_articles (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                project_id INTEGER NOT NULL,
                article_id INTEGER,
                external_id TEXT,
                title TEXT NOT NULL,
                authors TEXT,
                abstract TEXT,
                doi TEXT,
                pmid TEXT,
                url TEXT,
                publication_year INTEGER,
                journal TEXT,
                added_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                screening_status TEXT DEFAULT 'pending' CHECK (screening_status IN ('pending', 'included', 'excluded', 'conflict', 'resolved')),
                FOREIGN KEY (project_id) REFERENCES screening_projects(id) ON DELETE CASCADE,
                FOREIGN KEY (article_id) REFERENCES articles(id) ON DELETE SET NULL
            )
        """)

        # Screening Decisions Table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_decisions (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                screening_article_id INTEGER NOT NULL,
                reviewer_id INTEGER NOT NULL,
                decision TEXT NOT NULL CHECK (decision IN ('include', 'exclude', 'maybe')),
                reason TEXT,
                notes TEXT,
                screening_stage TEXT DEFAULT 'title_abstract' CHECK (screening_stage IN ('title_abstract', 'full_text')),
                decision_date DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (screening_article_id) REFERENCES screening_articles(id) ON DELETE CASCADE,
                FOREIGN KEY (reviewer_id) REFERENCES screening_reviewers(id) ON DELETE CASCADE,
                UNIQUE(screening_article_id, reviewer_id, screening_stage)
            )
        """)

        # Screening Conflicts Table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_conflicts (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                screening_article_id INTEGER NOT NULL,
                conflict_type TEXT NOT NULL CHECK (conflict_type IN ('decision_mismatch', 'quality_concern')),
                status TEXT DEFAULT 'unresolved' CHECK (status IN ('unresolved', 'resolved', 'escalated')),
                resolution_notes TEXT,
                resolved_by INTEGER,
                resolved_at DATETIME,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (screening_article_id) REFERENCES screening_articles(id) ON DELETE CASCADE,
                FOREIGN KEY (resolved_by) REFERENCES screening_reviewers(id) ON DELETE SET NULL
            )
        """)

        conn.commit()
        print("[DB] ✅ Database initialized successfully")


# -------------------------------------------------
# Connection Helper
# -------------------------------------------------
def get_connection():
    """Return a SQLite connection with row_factory set to sqlite3.Row."""
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
            cursor.execute("INSERT OR IGNORE INTO blog_tags (name) VALUES (?)", (tag,))
            cursor.execute("SELECT id FROM blog_tags WHERE name = ?", (tag,))
            tag_id = cursor.fetchone()[0]

            cursor.execute("""
                INSERT OR IGNORE INTO blog_tag_map (blog_id, tag_id)
                VALUES (?, ?)
            """, (blog_id, tag_id))

        conn.commit()
        print(f"[DB] ✅ Blog saved: {title} | Tags: {tags}")


def get_posts_by_tag(tag):
    """
    Fetch all blog posts that have the given tag.
    This uses proper JOINs with the tag tables instead of relying on LIKE.
    """
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT b.id, b.slug, b.title, b.author, b.content, b.created_at
            FROM blogs b
            JOIN blog_tag_map m ON b.id = m.blog_id
            JOIN blog_tags t ON m.tag_id = t.id
            WHERE LOWER(t.name) = LOWER(?)
            ORDER BY b.created_at DESC
        """, (tag,))
        return cursor.fetchall()


def get_blog_by_slug(slug):
    """Fetch a single blog post by slug, including its tags."""
    with get_connection() as conn:
        cursor = conn.cursor()

        # Fetch blog
        cursor.execute("""
            SELECT id, title, author, content, created_at
            FROM blogs
            WHERE slug = ?
        """, (slug,))
        row = cursor.fetchone()
        if not row:
            return None

        blog_id, title, author, content, created_at = row

        # Fetch tags
        cursor.execute("""
            SELECT name
            FROM blog_tags
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
    """Fetch all blog posts with their tags."""
    with get_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("""
            SELECT id, slug, title, author, created_at
            FROM blogs
            ORDER BY created_at DESC
        """)
        blogs = []
        for row in cursor.fetchall():
            blog_id, slug, title, author, created_at = row

            cursor.execute("""
                SELECT name
                FROM blog_tags
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
# Search & Article Storage
# -------------------------------------------------
def save_search(source, query, filters, results):
    """Save a search and its associated articles."""
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
    """Retrieve search history (latest first) with enhanced data structure."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT s.id, s.source, s.query, s.filters, s.timestamp,
                       COUNT(a.id) as article_count
                FROM searches s
                LEFT JOIN articles a ON s.id = a.search_id
                GROUP BY s.id, s.source, s.query, s.filters, s.timestamp
                ORDER BY s.timestamp DESC
            """)

            # Convert to dictionary format for easier handling
            rows = cursor.fetchall()
            history = []
            for row in rows:
                history.append({
                    'id': row[0],
                    'source': row[1],
                    'query': row[2],
                    'filters': row[3],
                    'timestamp': row[4],
                    'article_count': row[5],
                    'search_type': 'database'  # Default type for database searches
                })

            return history
    except Exception as e:
        print(f"[DB] ⚠ Failed to retrieve history: {e}")
        return []


def get_search_history_by_type(search_type='all', search_query=None, source_filter=None):
    """
    Retrieve search history categorized by type with filtering.

    Args:
        search_type (str): 'local', 'global', or 'all'
        search_query (str): Optional search term
        source_filter (str): Optional source filter

    Returns:
        list: Filtered and categorized search history
    """
    try:
        with get_connection() as conn:
            cursor = conn.cursor()

            # Base query with article count
            base_query = """
                SELECT s.id, s.source, s.query, s.filters, s.timestamp,
                       COUNT(a.id) as article_count
                FROM searches s
                LEFT JOIN articles a ON s.id = a.search_id
            """

            # Build WHERE conditions
            where_conditions = []
            params = []

            # Apply search filtering based on type
            if search_query:
                if search_type == 'local':
                    # Local search: only in query field
                    where_conditions.append("LOWER(s.query) LIKE ?")
                    params.append(f"%{search_query.lower()}%")
                elif search_type == 'global':
                    # Global search: across query, source, and timestamp
                    where_conditions.append("""
                        (LOWER(s.query) LIKE ? OR
                         LOWER(s.source) LIKE ? OR
                         LOWER(s.timestamp) LIKE ?)
                    """)
                    search_term = f"%{search_query.lower()}%"
                    params.extend([search_term, search_term, search_term])

            # Apply source filter
            if source_filter:
                where_conditions.append("s.source = ?")
                params.append(source_filter)

            # Construct final query
            where_clause = ""
            if where_conditions:
                where_clause = "WHERE " + " AND ".join(where_conditions)

            full_query = f"""
                {base_query}
                {where_clause}
                GROUP BY s.id, s.source, s.query, s.filters, s.timestamp
                ORDER BY s.timestamp DESC
            """

            cursor.execute(full_query, params)
            rows = cursor.fetchall()

            # Convert to enhanced dictionary format
            history = []
            for row in rows:
                history.append({
                    'id': row[0],
                    'source': row[1],
                    'query': row[2],
                    'filters': row[3],
                    'timestamp': row[4],
                    'article_count': row[5],
                    'search_type': search_type if search_type != 'all' else 'database',
                    'matched_fields': get_matched_fields(row, search_query, search_type) if search_query else []
                })

            return history

    except Exception as e:
        print(f"[DB] ⚠ Failed to retrieve search history by type: {e}")
        return []


def get_matched_fields(row, search_query, search_type):
    """
    Determine which fields matched the search query for highlighting.

    Args:
        row: Database row data
        search_query (str): The search term
        search_type (str): 'local' or 'global'

    Returns:
        list: Fields that matched the search
    """
    if not search_query:
        return []

    matched = []
    search_lower = search_query.lower()

    # Check query field (always checked for both local and global)
    if search_lower in str(row[2]).lower():
        matched.append('query')

    # For global search, also check source and timestamp
    if search_type == 'global':
        if search_lower in str(row[1]).lower():
            matched.append('source')
        if search_lower in str(row[4]).lower():
            matched.append('timestamp')

    return matched


def get_results_by_search_id(search_id):
    """Retrieve all articles from a given search ID."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path
                FROM articles
                WHERE search_id = ?
            """, (search_id,))
            return cursor.fetchall()
    except Exception as e:
        print(f"[DB] ⚠ Failed to retrieve results: {e}")
        return []


# Pagination and search functions removed as requested


# -------------------------------------------------
# Screening Project Functions
# -------------------------------------------------
def create_screening_project(name, description, inclusion_criteria, exclusion_criteria, created_by, require_dual_screening=False, conflict_resolution_method='discussion'):
    """Create a new screening project."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO screening_projects
                (name, description, inclusion_criteria, exclusion_criteria, created_by, require_dual_screening, conflict_resolution_method)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (name, description, inclusion_criteria, exclusion_criteria, created_by, require_dual_screening, conflict_resolution_method))
            project_id = cursor.lastrowid
            conn.commit()
            print(f"[DB] ✅ Screening project created: {name} (ID: {project_id})")
            return project_id
    except Exception as e:
        print(f"[DB] ⚠ Error creating screening project: {e}")
        return None


def get_screening_projects():
    """Get all screening projects."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT id, name, description, created_by, created_at, status, require_dual_screening
                FROM screening_projects
                ORDER BY created_at DESC
            """)
            return cursor.fetchall()
    except Exception as e:
        print(f"[DB] ⚠ Error retrieving screening projects: {e}")
        return []


def get_screening_project(project_id):
    """Get a specific screening project by ID."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM screening_projects WHERE id = ?
            """, (project_id,))
            return cursor.fetchone()
    except Exception as e:
        print(f"[DB] ⚠ Error retrieving screening project: {e}")
        return None


def add_reviewer_to_project(project_id, reviewer_name, reviewer_email, role='reviewer'):
    """Add a reviewer to a screening project."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO screening_reviewers (project_id, reviewer_name, reviewer_email, role)
                VALUES (?, ?, ?, ?)
            """, (project_id, reviewer_name, reviewer_email, role))
            reviewer_id = cursor.lastrowid
            conn.commit()
            print(f"[DB] ✅ Reviewer added to project {project_id}: {reviewer_name}")
            return reviewer_id
    except Exception as e:
        print(f"[DB] ⚠ Error adding reviewer: {e}")
        return None


def get_project_reviewers(project_id):
    """Get all reviewers for a project."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT id, reviewer_name, reviewer_email, role, added_at
                FROM screening_reviewers
                WHERE project_id = ?
                ORDER BY added_at
            """, (project_id,))
            return cursor.fetchall()
    except Exception as e:
        print(f"[DB] ⚠ Error retrieving project reviewers: {e}")
        return []


# -------------------------------------------------
# Screening Articles Functions
# -------------------------------------------------
def add_articles_to_screening_project(project_id, articles):
    """Add articles to a screening project."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()
            added_count = 0

            for article in articles:
                cursor.execute("""
                    INSERT INTO screening_articles
                    (project_id, article_id, external_id, title, authors, abstract, doi, pmid, url, publication_year, journal)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    project_id,
                    article.get('id'),
                    article.get('external_id'),
                    article.get('title', 'Untitled'),
                    article.get('authors', ''),
                    article.get('abstract', ''),
                    article.get('doi'),
                    article.get('pmid'),
                    article.get('url'),
                    article.get('publication_year'),
                    article.get('journal')
                ))
                added_count += 1

            conn.commit()
            print(f"[DB] ✅ Added {added_count} articles to screening project {project_id}")
            return added_count
    except Exception as e:
        print(f"[DB] ⚠ Error adding articles to screening project: {e}")
        return 0


def get_screening_articles(project_id, status=None, reviewer_id=None):
    """Get articles for screening in a project."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()

            base_query = """
                SELECT sa.id, sa.title, sa.authors, sa.abstract, sa.doi, sa.pmid,
                       sa.url, sa.publication_year, sa.journal, sa.screening_status,
                       sa.added_at
                FROM screening_articles sa
                WHERE sa.project_id = ?
            """
            params = [project_id]

            if status:
                base_query += " AND sa.screening_status = ?"
                params.append(status)

            if reviewer_id:
                base_query += """
                    AND sa.id NOT IN (
                        SELECT sd.screening_article_id
                        FROM screening_decisions sd
                        WHERE sd.reviewer_id = ?
                    )
                """
                params.append(reviewer_id)

            base_query += " ORDER BY sa.added_at"

            cursor.execute(base_query, params)
            return cursor.fetchall()
    except Exception as e:
        print(f"[DB] ⚠ Error retrieving screening articles: {e}")
        return []


def save_screening_decision(screening_article_id, reviewer_id, decision, reason=None, notes=None, screening_stage='title_abstract'):
    """Save a screening decision."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()

            # Insert or update decision
            cursor.execute("""
                INSERT OR REPLACE INTO screening_decisions
                (screening_article_id, reviewer_id, decision, reason, notes, screening_stage)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (screening_article_id, reviewer_id, decision, reason, notes, screening_stage))

            # Update article screening status
            update_article_screening_status(cursor, screening_article_id)

            conn.commit()
            print(f"[DB] ✅ Screening decision saved: {decision} for article {screening_article_id}")
            return True
    except Exception as e:
        print(f"[DB] ⚠ Error saving screening decision: {e}")
        return False


def update_article_screening_status(cursor, screening_article_id):
    """Update the overall screening status of an article based on decisions."""
    # Get all decisions for this article
    cursor.execute("""
        SELECT decision FROM screening_decisions
        WHERE screening_article_id = ?
    """, (screening_article_id,))
    decisions = [row[0] for row in cursor.fetchall()]

    if not decisions:
        status = 'pending'
    elif len(decisions) == 1:
        status = 'included' if decisions[0] == 'include' else 'excluded'
    else:
        # Multiple decisions - check for conflicts
        unique_decisions = set(decisions)
        if len(unique_decisions) == 1:
            status = 'included' if decisions[0] == 'include' else 'excluded'
        else:
            status = 'conflict'

    cursor.execute("""
        UPDATE screening_articles
        SET screening_status = ?
        WHERE id = ?
    """, (status, screening_article_id))


def get_screening_statistics(project_id):
    """Get screening statistics for a project."""
    try:
        with get_connection() as conn:
            cursor = conn.cursor()

            # Total articles
            cursor.execute("SELECT COUNT(*) FROM screening_articles WHERE project_id = ?", (project_id,))
            total_articles = cursor.fetchone()[0]

            # Articles by status
            cursor.execute("""
                SELECT screening_status, COUNT(*)
                FROM screening_articles
                WHERE project_id = ?
                GROUP BY screening_status
            """, (project_id,))
            status_counts = dict(cursor.fetchall())

            # Reviewer progress
            cursor.execute("""
                SELECT sr.reviewer_name, COUNT(sd.id) as decisions_made
                FROM screening_reviewers sr
                LEFT JOIN screening_decisions sd ON sr.id = sd.reviewer_id
                WHERE sr.project_id = ?
                GROUP BY sr.id, sr.reviewer_name
            """, (project_id,))
            reviewer_progress = cursor.fetchall()

            return {
                'total_articles': total_articles,
                'status_counts': status_counts,
                'reviewer_progress': reviewer_progress
            }
    except Exception as e:
        print(f"[DB] ⚠ Error getting screening statistics: {e}")
        return None

