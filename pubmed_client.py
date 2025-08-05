# pubmed_client.py

from Bio import Entrez
import time
#import requests

# Set Entrez credentials
Entrez.email = "Your_Email.com"
Entrez.api_key = "Your_API_key"


# ✅ Throttle: NCBI recommends <= 3 requests/sec (0.34s delay for safety)
THROTTLE_DELAY = 0.34


def search_pubmed(query, author=None, start_year=None, end_year=None, country=None, max_results=10):
    """
    PubMed: Search articles with keyword/boolean/DOI + filters (author/year/country).
    Returns full details including DOI, PMID, Abstract, Authors, Affiliations.
    """

    # ✅ Build search query
    search_query = query.strip()

    # Author filter
    if author:
        search_query += f' AND {author}[Author]'

    # Year range filter
    if start_year or end_year:
        start = start_year or "1800"
        end = end_year or "3000"
        search_query += f' AND ("{start}"[Date - Publication] : "{end}"[Date - Publication])'

    # Country filter (via affiliations)
    if country:
        search_query += f' AND {country}[Affiliation]'

    # ✅ Search PubMed IDs
    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record.get("IdList", [])
    time.sleep(THROTTLE_DELAY)

    if not ids:
        return []

    # ✅ Fetch details for IDs
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="xml")
    papers = Entrez.read(handle)
    time.sleep(THROTTLE_DELAY)

    results = []
    for paper in papers.get("PubmedArticle", []):
        article_data = paper.get("MedlineCitation", {}).get("Article", {})

        # Title
        title = article_data.get("ArticleTitle", "No title")

        # Authors
        authors_data = article_data.get("AuthorList", [])
        authors = [
            f"{a.get('ForeName', '')} {a.get('LastName', '')}".strip()
            for a in authors_data if a.get("ForeName") or a.get("LastName")
        ]

        # Affiliations
        affiliations = []
        for author_entry in authors_data:
            for aff in author_entry.get("AffiliationInfo", []):
                aff_text = aff.get("Affiliation", "").strip()
                if aff_text:
                    affiliations.append(aff_text)

        # DOI
        doi = None
        article_ids = paper.get("PubmedData", {}).get("ArticleIdList", [])
        for article_id in article_ids:
            if article_id.attributes.get("IdType") == "doi":
                doi = str(article_id)

        # Abstract
        abstract_list = article_data.get("Abstract", {}).get("AbstractText", [])
        abstract_text = " ".join(abstract_list) if abstract_list else "No abstract available"

        # PMID & URL
        pmid = paper.get("MedlineCitation", {}).get("PMID", "")
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        # ✅ Append result in DB-compatible format
        results.append({
            "title": title,
            "authors": authors,
            "abstract": abstract_text,
            "affiliations": affiliations,
            "pmid": pmid,
            "doi": doi,
            "url": url,
            "source": "PubMed"
        })

    return results

