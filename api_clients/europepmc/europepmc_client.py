#europepmc

# europepmc_client.py
import requests

BASE_URL_REST = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
BASE_URL_GRIST = "https://www.ebi.ac.uk/europepmc/GristAPI/rest/get"

def search_epmc_publications(query, page=1, format_type="json"):
    """Search Europe PMC REST API for publications (returns parsed dict list)."""
    params = {
        "query": query,
        "format": format_type,
        "page": page,
        "resultType": "core"  # Core includes abstracts, affiliations, etc.
    }
    
    response = requests.get(BASE_URL_REST, params=params)
    
    if response.status_code != 200:
        return [{"error": f"Request failed with status {response.status_code}"}]
    
    data = response.json()
    results = data.get("resultList", {}).get("result", [])
    
    parsed_results = []
    for item in results:
        parsed_results.append({
            "title": item.get("title", "Not available"),
            "authors": item.get("authorString", "").split(", ") if item.get("authorString") else [],
            "doi": item.get("doi", "Not available"),
            "abstract": item.get("abstractText", "Not available"),
            "affiliations": [item.get("journalTitle", "Not available")],
            "pmid": item.get("pmid", None),
            "url": f"https://europepmc.org/article/MED/{item.get('pmid')}" if item.get("pmid") else "Not available"
        })
    
    return parsed_results

def search_epmc_grants(query, page=1, result_type="core", format_type="json"):
    """Search Europe PMC Grist API for grants (returns parsed dict list)."""
    params = {
        "query": query,
        "page": page,
        "resultType": result_type,
        "format": format_type
    }
    
    headers = {"Accept": "application/json"}
    response = requests.get(BASE_URL_GRIST, params=params, headers=headers)
    
    if response.status_code != 200:
        return [{"error": f"Request failed with status {response.status_code}"}]
    
    data = response.json()
    results = data.get("resultList", {}).get("result", [])
    
    parsed_results = []
    for item in results:
        parsed_results.append({
            "title": item.get("title", "Not available"),
            "authors": [item.get("pi", "Not available")],
            "doi": None,  # Grants usually don’t have DOI
            "abstract": item.get("abstract", "Not available"),
            "affiliations": [item.get("aff", "Not available")],
            "pmid": None,
            "url": "Not available"
        })
    
    return parsed_results

