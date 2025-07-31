
# corssref_crossref_client.py

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_URL = "https://api.crossref.org/works"

# 🔄 Session with Retry for better reliability
session = requests.Session()
retries = Retry(total=3, backoff_factor=2, status_forcelist=[500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retries))

def search_crossref(query, rows=10):
    """Crossref: Keyword search for minimal metadata."""
    params = {"query": query, "rows": rows}
    try:
        r = session.get(BASE_URL, params=params, timeout=20)  # ⏳ Increased timeout
        r.raise_for_status()
    except requests.exceptions.Timeout:
        print("⚠ Crossref API timed out. Try a more specific keyword.")
        return []
    except requests.exceptions.RequestException as e:
        print(f"⚠ Crossref API error: {e}")
        return []

    items = r.json().get("message", {}).get("items", [])
    results = []
    for item in items:
        title = item.get("title", ["No title available"])[0]
        authors = [
            f"{a.get('given', '')} {a.get('family', '')}".strip()
            for a in item.get("author", [])
        ]
        doi = item.get("DOI", "")
        url = item.get("URL", "")
        results.append({
            "title": title,
            "authors": authors,
            "doi": doi,
            "url": url,
            "source": "Crossref"
        })
    return results

