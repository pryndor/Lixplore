# springer_client
import requests
import os

# 🔑 Springer API key (can also set via environment variable)
API_KEY = os.getenv("SPRINGER_API_KEY", "2849e4acd5190c3992f9a1a8ebe82692")

# 📂 Folder to save PDFs
SAVE_FOLDER = "springer_pdfs"
os.makedirs(SAVE_FOLDER, exist_ok=True)


def search_springer(query, results=5, author=None, start_year=None, end_year=None):
    """
    Springer: Search articles by keyword or DOI.
    Returns metadata + abstract + OA PDF link (if available).
    """

    # ✅ Add year filter
    if start_year and end_year:
        query += f" AND date-facet-start:[{start_year} TO {end_year}]"
    elif start_year:
        query += f" AND date-facet-start:[{start_year} TO 3000]"
    elif end_year:
        query += f" AND date-facet-start:[1800 TO {end_year}]"

    metadata_url = f"https://api.springernature.com/metadata/json?q={query}&p={results}&api_key={API_KEY}"
    print(f"\n🔍 Searching Springer for: {query}")

    try:
        response = requests.get(metadata_url, timeout=15)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"❌ Springer API error: {e}")
        return []

    data = response.json()
    if not data.get("records"):
        print("⚠ No results found.")
        return []

    results_list = []
    for record in data["records"]:
        title = record.get("title", "No title")
        abstract = record.get("abstract", "No abstract available")
        doi = record.get("doi")
        authors = [a.get("creator") for a in record.get("creators", []) if a.get("creator")]

        # ✅ URL: Prefer HTML link, else construct from DOI
        url_value = None
        for url_item in record.get("url", []):
            if url_item.get("format") == "html":
                url_value = url_item.get("value")
                break
        if not url_value and doi:
            url_value = f"https://doi.org/{doi}"

        # ✅ OA PDF (optional download)
        oa_pdf_path = None
        for url_item in record.get("url", []):
            if url_item.get("format") == "pdf":
                oa_pdf_path = url_item.get("value")
                download_pdf(oa_pdf_path, f"{SAVE_FOLDER}/{title[:50].replace('/', '_')}.pdf")
                break

        results_list.append({
            "title": title,
            "authors": authors,
            "abstract": abstract,
            "affiliations": [],
            "pmid": None,
            "doi": doi,
            "url": url_value,
            "oa_pdf_path": oa_pdf_path,
            "source": "Springer"
        })

    return results_list

