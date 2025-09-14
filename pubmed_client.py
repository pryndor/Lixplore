# pubmed_client.py

import requests
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from dotenv import load_dotenv
import os

# Load .env variables
load_dotenv()

# Environment variables (ensure these exist in your .env file)
NCBI_EMAIL = os.getenv("NCBI_EMAIL")    # .env: NCBI_EMAIL=YOUR_Email
NCBI_API_KEY = os.getenv("NCBI_API_KEY")   # .env: NCBI_API_KEY=YOUR_API
PMC_CONVERTER_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"

# API base URLs
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"




# Throttle delay per NCBI guidelines (max 3 requests/sec)
THROTTLE_DELAY = 0.125


# ------------------------
# 1) Old Entrez-based search
# ------------------------
def search_pubmed_entrez(query):
    """
    Search PubMed using Biopython's Entrez API.
    Returns a list of PMIDs.
    """
    from Bio import Entrez
    Entrez.email = NCBI_EMAIL
    if NCBI_API_KEY:
        Entrez.api_key = NCBI_API_KEY

    handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
    record = Entrez.read(handle)
    return record["IdList"]


# ------------------------
# 2) Direct HTTP-based search with filters
# ------------------------
def search_pubmed(query, author=None, start_year=None, end_year=None, country=None, max_results=20):
    """
    Search PubMed using direct HTTP requests to NCBI E-utilities.
    Filters: author, publication year range, country (affiliation).
    Returns list of article dicts with metadata, including PMC info if available.
    """

    # Build query string
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

    # 1) ESearch: get PMIDs
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

    root = ET.fromstring(esearch_resp.text)
    idlist = root.find("IdList")
    if idlist is None:
        return []
    pmids = [id_elem.text for id_elem in idlist.findall("Id")]
    if not pmids:
        return []

    # 2) EFetch: fetch article details
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

    root = ET.fromstring(efetch_resp.text)
    articles = []

    for article_elem in root.findall(".//PubmedArticle"):
        medline_citation = article_elem.find("MedlineCitation")
        article = medline_citation.find("Article") if medline_citation is not None else None
        if article is None:
            continue

        # Title
        title = article.findtext("ArticleTitle", default="No title")

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
        abstract = article.find("Abstract")
        if abstract is not None:
            abstract_text = " ".join([elem.text or "" for elem in abstract.findall("AbstractText")]).strip()
        else:
            abstract_text = "No abstract available"

        # PMID
        pmid = medline_citation.findtext("PMID", default="")

        # PubMed URL
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None

        # PMC info
        pmc_result = get_pmc_pdf_url(pmid)

        articles.append({
            "title": title,
            "authors": authors,
            "abstract": abstract_text,
            "affiliations": affiliations,
            "pmid": pmid,
            "doi": doi,
            "url": url,
            "source": "PubMed",
            "pmcid": pmc_result["pmcid"] if pmc_result else None,
            "pmc_pdf_url": pmc_result["pdf_url"] if pmc_result else None
        })

    return articles


# ---------- PMID -> PMCID PDF helper ----------


def get_pmc_pdf_url(pmid):
    """
    Given a PubMed ID, returns the PMC ID and direct PDF URL if available.
    """
    params = {
        "ids": pmid,
        "tool": "Lixplore",
        "email": NCBI_EMAIL,
        "format": "json"
    }

    response = requests.get(PMC_CONVERTER_URL, params=params)
    time.sleep(THROTTLE_DELAY)

    if response.status_code != 200:
        return None

    data = response.json()
    if "records" not in data or not data["records"]:
        return None

    record = data["records"][0]
    pmcid = record.get("pmcid")
    if not pmcid:
        return None

    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/"
    return {"pmcid": pmcid, "pdf_url": pdf_url}

