# pubmed_client.py

from Bio import Entrez
import time

# Set Entrez credentials
Entrez.email = "balathepharmacist@gmail.com"
Entrez.api_key = "781f12bc04105f7d5536a510520cd74cbf08"

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






