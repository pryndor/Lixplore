from Bio import Entrez
import time

Entrez.email = "balathepharmacist@gmail.com"
Entrez.api_key = "781f12bc04105f7d5536a510520cd74cbf08"

def search_pubmed(query, max_results=5):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]

    results = []

    for pmid in ids:
        time.sleep(0.1)  # limit to 10 req/sec
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
        fetch_record = Entrez.read(fetch_handle)
        fetch_handle.close()

        for article in fetch_record["PubmedArticle"]:
            article_data = article["MedlineCitation"]["Article"]
            title = article_data.get("ArticleTitle", "No title available")
            abstract = article_data.get("Abstract", {}).get("AbstractText", ["No abstract available"])[0]

            results.append({
                "title": title,
                "abstract": abstract
            })

    return results

