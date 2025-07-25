# pubmed_client.py

from Bio import Entrez
import time
import requests

# Set Entrez credentials
Entrez.email = "Your_email.com"
Entrez.api_key = "Your_API_Key"

def search_pubmed(query, max_results=5, mindate=None, maxdate=None, country=None):
    # Prepare search arguments
    search_args = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "xml"
    }

    if mindate:
        search_args["mindate"] = mindate
    if maxdate:
        search_args["maxdate"] = maxdate

    # Execute search
    handle = Entrez.esearch(**search_args)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])

    results = []

    for pmid in ids:
        time.sleep(0.1)  # Respect NCBI rate limits
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
        fetch_record = Entrez.read(fetch_handle)
        fetch_handle.close()

        for article in fetch_record.get("PubmedArticle", []):
            citation = article.get("MedlineCitation", {})
            article_data = citation.get("Article", {})

            title = article_data.get("ArticleTitle", "No title available")

            # Handle abstract as string (join parts if needed)
            abstract_parts = article_data.get("Abstract", {}).get("AbstractText", [])
            if isinstance(abstract_parts, list):
                abstract = " ".join(abstract_parts)
            else:
                abstract = abstract_parts or "No abstract available"

            # Extract affiliations
            affiliations = []
            authors = article_data.get("AuthorList", [])
            for author in authors:
                affs = author.get("AffiliationInfo", [])
                for aff in affs:
                    aff_text = aff.get("Affiliation")
                    if aff_text:
                        affiliations.append(aff_text)

            # Filter by country if specified
            if country:
                country_lower = country.lower().strip()
                matched = any(country_lower in (aff or "").lower() for aff in affiliations)
                if not matched:
                    continue  # Skip this article

            # Add result
            results.append({
                "title": title,
                "abstract": abstract,
                "affiliations": affiliations
            })

    return results

import requests

def fetch_pmc_pdf(pmid: str):
    """
    Given a PMID, fetch and download the PMC full-text PDF if available.
    """
    try:
        link_handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        link_result = Entrez.read(link_handle)
        link_handle.close()

        linksets = link_result[0].get("LinkSetDb")
        if not linksets:
            print("⚠️ No PMC full-text found for this PMID.")
            return

        pmc_id = linksets[0]["Link"][0]["Id"]
        pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/"
        filename = f"PMC{pmc_id}.pdf"

        r = requests.get(pdf_url, stream=True, timeout=10)
        if r.status_code == 200 and "application/pdf" in r.headers.get("Content-Type", ""):
            with open(filename, "wb") as f:
                for chunk in r.iter_content(1024):
                    f.write(chunk)
            print(f"✅ PDF downloaded: {filename}")
        else:
            print(f"❌ PDF not available. HTTP status: {r.status_code}")
    except Exception as e:
        print(f"❌ Error fetching PDF: {e}")






