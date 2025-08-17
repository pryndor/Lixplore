# springer_client

# springer_client_refactored.py
from springernature_api_client.openaccess import OpenAccessAPI
from springernature_api_client.metadata import MetadataAPI
from springernature_api_client.tdm import TDMAPI
from springernature_api_client.utils import results_to_dataframe
import os
import requests

# 🔑 API key
API_KEY = os.getenv("SPRINGER_API_KEY", "2849e4acd5190c3992f9a1a8ebe82692")

# 📂 Folder to save PDFs
SAVE_FOLDER = "springer_pdfs"
os.makedirs(SAVE_FOLDER, exist_ok=True)


def download_pdf(url, filepath):
    """Download PDF if URL exists"""
    if url:
        try:
            r = requests.get(url)
            if r.status_code == 200:
                with open(filepath, "wb") as f:
                    f.write(r.content)
                print(f"✅ PDF saved: {filepath}")
            else:
                print(f"❌ Failed to download PDF: {url}")
        except Exception as e:
            print(f"❌ Error downloading PDF: {e}")


def search_springer(query, results=5):
    """
    Search Springer articles using Open Access API, Metadata API, and TDM API.
    Returns list of dicts with metadata + abstract + OA PDF link.
    """
    oa_client = OpenAccessAPI(api_key=API_KEY)
    metadata_client = MetadataAPI(api_key=API_KEY)
    tdm_client = TDMAPI(api_key=API_KEY)  # optional, for full-text access

    print(f"\n🔍 Searching Springer for: {query}")
    response = oa_client.search(q=f'keyword:"{query}"', p=results, s=1, fetch_all=False, is_premium=False)

    articles = []
    for record in response.get("records", []):
        title = record.get("title", "No title")
        abstract = record.get("abstract", "No abstract")
        doi = record.get("doi")
        authors = [a.get("creator") for a in record.get("creators", []) if a.get("creator")]

        # HTML URL fallback
        url_value = None
        for url_item in record.get("url", []):
            if url_item.get("format") == "html":
                url_value = url_item.get("value")
                break
        if not url_value and doi:
            url_value = f"https://doi.org/{doi}"

        # OA PDF
        oa_pdf_path = None
        for url_item in record.get("url", []):
            if url_item.get("format") == "pdf":
                oa_pdf_path = url_item.get("value")
                download_pdf(oa_pdf_path, f"{SAVE_FOLDER}/{title[:50].replace('/', '_')}.pdf")
                break

        articles.append({
            "title": title,
            "authors": authors,
            "abstract": abstract,
            "doi": doi,
            "url": url_value,
            "oa_pdf_path": oa_pdf_path,
            "source": "Springer"
        })

    return articles

if __name__ == "__main__":
    # Test fetching PDFs
    test_articles = search_springer("cancer", results=5)
    for article in test_articles:
        if article["oa_pdf_path"]:
            print("OA PDF available:", article["oa_pdf_path"])
            # Optional: download
            # download_pdf(article["oa_pdf_path"], f"{SAVE_FOLDER}/{article['title'][:50].replace('/', '_')}.pdf")
        else:
            print("No PDF available for:", article["title"])


