# main.py


from pubmed_client import search_pubmed
from crossref.crossref_client import get_sample_crossref_result

def main():
    print("🔍 PubMed Literature Search Tool\n")

    # Choose search type
    search_type = input("Search by: [1] Keyword  [2] DOI  [3] PMID : ").strip()

    if search_type == "2":
        doi = input("Enter DOI: ").strip()
        query = f"{doi}[DOI]"
    elif search_type == "3":
        pmid = input("Enter PMID: ").strip()
        query = f"{pmid}[PMID]"
    else:
        # Default to keyword search
        keywords = input("Enter keywords (comma-separated): ").strip()
        operator = input("Choose operator [AND/OR/NOT] (optional, press Enter to skip): ").strip().upper()
        if operator in ["AND", "OR", "NOT"]:
            query = f" {operator} ".join([k.strip() for k in keywords.split(",")])
        else:
            query = " ".join([k.strip() for k in keywords.split(",")])

    # Optional filters
    from_year = input("From year (optional): ").strip()
    to_year = input("To year (optional): ").strip()
    country = input("Country filter (optional): ").strip()

    print(f"\n🔍 Searching PubMed for: {query}\n")

    # Fetch more than needed to allow rejections (say 20)
    results = search_pubmed(
        query,
        max_results=20,
        mindate=from_year if from_year else None,
        maxdate=to_year if to_year else None,
        country=country if country else None
    )

    if not results:
        print("❌ No results found.")
        return

    shown = 0
    index = 0
    accepted = 0
    max_show = 10

    while accepted < max_show and index < len(results):
        article = results[index]
        print(f"\n📝 Result {index + 1}:")
        print(f"📌 Title: {article['title']}")
        print(f"📄 Abstract: {article['abstract']}")
        if article.get("affiliations"):
            print("🏥 Affiliations:")
            for aff in article["affiliations"]:
                print(f"  - {aff}")

        choice = input("\nDo you want to keep this article? [Y/n]: ").strip().lower()
        if choice in ["", "y", "yes"]:
            accepted += 1
        else:
            print("⏭️ Skipped.")

        index += 1

    if accepted == 0:
        print("\n😕 You skipped all articles.")
    else:
        print(f"\n✅ You accepted {accepted} articles.")

if __name__ == "__main__":
    main()

