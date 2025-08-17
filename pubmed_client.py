# pubmed_client.py

import requests
import time
import xml.etree.ElementTree as ET
from datetime import datetime
import os
from dotenv import load_dotenv

load_dotenv()

# NCBI API base URLs
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

PUBMED_EMAIL = os.getenv("balathepharmacist@gmail.com")
PUBMED_API_KEY = os.getenv("781f12bc04105f7d5536a510520cd74cbf08")


# Throttle delay per NCBI guidelines (max 3 requests/sec)
THROTTLE_DELAY = 0.34


def search_pubmed(query, author=None, start_year=None, end_year=None, country=None, max_results=10):
    """
    Search PubMed using direct HTTP requests to NCBI E-utilities.
    Filters: author, publication year range, country (affiliation).
    Returns list of article dicts with metadata.
    """

    # Build base query string with filters
    search_terms = [query.strip()]

    if author:
        search_terms.append(f'{author}[Author]')

    if start_year or end_year:
        start = start_year or "1800"
        end = end_year or datetime.now().year
        search_terms.append(f'("{start}"[Date - Publication] : "{end}"[Date - Publication])')

    if country:
        search_terms.append(f'{country}[Affiliation]')

    full_query = " AND ".join(search_terms)

    # 1) Use esearch to get list of PMIDs matching query
    esearch_params = {
        "db": "pubmed",
        "term": full_query,
        "retmax": max_results,
        "retmode": "xml",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    esearch_resp = requests.get(ESEARCH_URL, params=esearch_params)
    time.sleep(THROTTLE_DELAY)
    if esearch_resp.status_code != 200:
        print(f"Error in esearch request: {esearch_resp.status_code}")
        return []

    # Parse esearch XML response to get list of IDs
    root = ET.fromstring(esearch_resp.text)
    idlist = root.find("IdList")
    if idlist is None:
        return []
    pmids = [id_elem.text for id_elem in idlist.findall("Id")]
    if not pmids:
        return []

    # 2) Use efetch to get detailed info for each PMID
    efetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }
    efetch_resp = requests.get(EFETCH_URL, params=efetch_params)
    time.sleep(THROTTLE_DELAY)
    if efetch_resp.status_code != 200:
        print(f"Error in efetch request: {efetch_resp.status_code}")
        return []

    # Parse efetch XML for article details
    root = ET.fromstring(efetch_resp.text)
    articles = []

    for article_elem in root.findall(".//PubmedArticle"):
        medline_citation = article_elem.find("MedlineCitation")
        article = medline_citation.find("Article") if medline_citation is not None else None
        if article is None:
            continue

        # Title
        title_elem = article.find("ArticleTitle")
        title = title_elem.text if title_elem is not None else "No title"

        # Authors
        authors = []
        author_list = article.find("AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                fore_name = author.findtext("ForeName") or ""
                last_name = author.findtext("LastName") or ""
                full_name = (fore_name + " " + last_name).strip()
                if full_name:
                    authors.append(full_name)

        # Affiliations
        affiliations = []
        if author_list is not None:
            for author in author_list.findall("Author"):
                for aff_info in author.findall("AffiliationInfo"):
                    aff_text = aff_info.findtext("Affiliation")
                    if aff_text and aff_text not in affiliations:
                        affiliations.append(aff_text.strip())

        # DOI
        doi = None
        pubmed_data = article_elem.find("PubmedData")
        if pubmed_data is not None:
            article_ids = pubmed_data.find("ArticleIdList")
            if article_ids is not None:
                for article_id in article_ids.findall("ArticleId"):
                    if article_id.attrib.get("IdType") == "doi":
                        doi = article_id.text
                        break

        # Abstract
        abstract_text = ""
        abstract = article.find("Abstract")
        if abstract is not None:
            abstract_text = " ".join([elem.text or "" for elem in abstract.findall("AbstractText")]).strip()
        if not abstract_text:
            abstract_text = "No abstract available"

        # PMID
        pmid_elem = medline_citation.find("PMID") if medline_citation is not None else None
        pmid = pmid_elem.text if pmid_elem is not None else ""

        # URL to PubMed article
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None

        articles.append({
            "title": title,
            "authors": authors,
            "abstract": abstract_text,
            "affiliations": affiliations,
            "pmid": pmid,
            "doi": doi,
            "url": url,
            "source": "PubMed"
        })

    return articles

