# main.py

from pubmed_client import search_pubmed 
from crossref.crossref_client import get_sample_crossref_result

def main():
    print("🔍 PubMed Literature Search Tool")
    query = input("Enter search term: ")
    print(f"\nSearching for: {query}...\n")

    results = search_pubmed(query)

    if not results:
        print("❌ No results found.")
        return

    for i, r in enumerate(results, 1):
        print(f"🔹 Result {i}")
        print(f"Title   : {r['title']}")
        print(f"Abstract: {r['abstract']}\n")

    # Crossref (mock)
    crossref_result = get_sample_crossref_result()
    print("\n📘 Crossref Result:")
    print(f"\nTitle: {crossref_result['title']}")
    print(f"Abstract: {crossref_result['abstract']}")
    print(f"Authors: {', '.join(crossref_result['authors'])}")
    print(f"DOI: {crossref_result['doi']}")
    print(f"Source: {crossref_result['source']}")

if __name__ == "__main__":
    main()

