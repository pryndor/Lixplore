# DOAJ_client

import requests

BASE_URL = "https://doaj.org/api/v4/search/articles/"

def search_doaj(query, page=1, page_size=10):
    """
    Search DOAJ articles and return detailed information.
    
    :param query: Search keyword(s)
    :param page: Page number
    :param page_size: Results per page
    :return: List of article dicts
    """
    url = f"{BASE_URL}{query}"
    params = {
        "page": page,
        "pageSize": page_size
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        articles = []
        for result in data.get("results", []):
            bibjson = result.get("bibjson", {})
            
            # Title
            title = bibjson.get("title", "No title available")
            
            # Authors
            authors = [a.get("name") for a in bibjson.get("author", [])] or ["Author(s) not available"]
            
            # Journal
            journal = bibjson.get("journal", {}).get("title", "Unknown Journal")
            
            # DOI
            doi = next(
                (identifier.get("id") for identifier in bibjson.get("identifier", [])
                 if identifier.get("type", "").lower() == "doi"),
                None
            )
            
            # Abstract
            abstract = bibjson.get("abstract", "Abstract not available")
            
            # Article link
            article_url = next(
                (link.get("url") for link in bibjson.get("link", []) if link.get("url")),
                "No URL available"
            )
            
            articles.append({
                "title": title,
                "authors": authors,
                "journal": journal,
                "doi": doi or "DOI not available",
                "abstract": abstract,
                "url": article_url
            })
        
        return articles
    
    except requests.exceptions.RequestException as e:
        print(f"❌ DOAJ API request error: {e}")
        return []

