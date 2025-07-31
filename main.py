# main.py

from pubmed_client import search_pubmed
# from pdf_fetcher import fetch_pmc_pdf
from crossref.crossref_client import search_crossref
from database.database import init_db, save_search, get_archive
#from database.database import init_db, save_search  # ✅ Added database imports
# from scholarly_client import search_scholar  # Optional for later
from springer.springer_client import search_springer   # ✅ Added Springer import
from database.database import init_db, save_search, get_archive

from pubmed_client import search_pubmed
from crossref.crossref_client import search_crossref
from springer.springer_client import search_springer
from europepmc.europepmc_client import search_epmc_publications, search_epmc_grants
from database.database import init_db, save_search, get_archive
import json


def main():
    # ✅ Initialize database
    init_db()

    while True:
        print("\n🔍 Literature Search Tool\n")
        print("Select Option:")
        print("1. Crossref (Keyword search only: Metadata)")
        print("2. PubMed (Keyword / Boolean / Author / DOI search: Full data)")
        print("3. Springer (Keyword / Boolean / DOI search: Metadata + Abstract + OA PDFs)")
        print("4. Europe PMC (Publications & Grants)")
        print("5. View Search History")
        print("6. Exit")
        choice = input("Enter choice (1-6): ").strip()

        if choice == "6":
            print("\n👋 Exiting tool. Goodbye!")
            break

        filters = {}
        results = []
        source_name = ""

        # ✅ View Search History
        if choice == "5":
            history = get_archive()
            if not history:
                print("\n⚠ No search history found.")
            else:
                print("\n📜 Search History:")
                for h in history:
                    search_id, source, query, filt, timestamp = h
                    try:
                        filt_dict = eval(filt) if isinstance(filt, str) else filt
                        filt_text = ", ".join(f"{k}: {v}" for k, v in filt_dict.items() if v)
                    except Exception:
                        filt_text = str(filt)
                    print(f"[{timestamp}] Source: {source} | Query: {query} | Filters: {filt_text}")
            continue

        # ✅ PubMed Search
        elif choice == "2":
            print("\nSelect PubMed Search Type:")
            print("1. Keyword search")
            print("2. Boolean search (e.g., brain AND cancer, AI OR bias)")
            print("3. DOI search")
            print("4. Author search")
            search_type = input("Enter choice (1-4): ").strip()

            if search_type == "3":
                query = input("Enter DOI: ").strip()
            elif search_type == "4":
                author = input("Enter author name: ").strip()
                query = f"{author}[Author]"
                filters["author"] = author
            else:
                query = input("Enter keyword or Boolean query: ").strip()
                author = input("Author name (optional): ").strip()
                if author:
                    query += f" AND {author}[Author]"
                    filters["author"] = author

            start_year = input("Start year (optional): ").strip() or None
            end_year = input("End year (optional): ").strip() or None
            country = input("Country (optional): ").strip() or None
            filters.update({"start_year": start_year, "end_year": end_year, "country": country})

            results = search_pubmed(query, None, start_year, end_year, country)
            source_name = "PubMed"

        # ✅ Crossref Search
        elif choice == "1":
            query = input("Enter keyword: ").strip()
            results = search_crossref(query)
            results = [r for r in results if isinstance(r, dict)]
            filters["note"] = "Keyword search only (Crossref metadata)"
            source_name = "Crossref"

        # ✅ Springer Search
        elif choice == "3":
            print("\nSelect Springer Search Type:")
            print("1. Keyword search")
            print("2. Boolean search (use AND / OR)")
            print("3. DOI search")
            search_type = input("Enter choice (1-3): ").strip()

            if search_type == "3":
                query = input("Enter DOI: ").strip()
            elif search_type == "2":
                print("\n⚠ Example: pharmacovigilance AND safety OR risk")
                query = input("Enter Boolean query: ").strip()
            else:
                query = input("Enter keyword: ").strip()

            start_year = input("Start year (optional): ").strip() or None
            end_year = input("End year (optional): ").strip() or None
            if start_year or end_year:
                filters.update({"start_year": start_year, "end_year": end_year})

            results = search_springer(query, results=5, start_year=start_year, end_year=end_year)
            source_name = "Springer"

        # ✅ Europe PMC Search
        elif choice == "4":
            print("\nSelect Europe PMC Search Type:")
            print("1. Publications (Articles, Books, Preprints)")
            print("2. Grants (Funding Information)")
            epmc_type = input("Enter choice (1-2): ").strip()

            query = input("Enter keyword or fielded query: ").strip()
            page = input("Page number (optional, default=1): ").strip() or 1

            if epmc_type == "1":
                results = search_epmc_publications(query, page=page)
                source_name = "Europe PMC Publications"
            elif epmc_type == "2":
                results = search_epmc_grants(query, page=page)
                source_name = "Europe PMC Grants"
            else:
                print("\n⚠ Invalid choice for Europe PMC.")
                continue

        else:
            print("\n⚠ Invalid choice. Try again.")
            continue

        if not results:
            print("\n⚠ No results found.")
            continue

        # ✅ Save search & results
        try:
            save_search(source_name, query, json.dumps(filters), json.dumps(results))
        except Exception as e:
            print(f"⚠ Failed to save search: {e}")

        # ✅ Display results
        print(f"\n✅ Found {len(results)} results:")
        for idx, article in enumerate(results, start=1):
            print(f"\nResult {idx} ({source_name})")
            print(f"Title: {article.get('title', 'Not available')}")
            print(f"Authors: {', '.join(article.get('authors', [])) if article.get('authors') else 'Not available'}")
            print(f"DOI: {article.get('doi', 'Not available')}")
            if source_name == "PubMed":
                print(f"PMID: {article.get('pmid', 'Not available')}")
            print(f"Abstract: {article.get('abstract', 'Not available')}")
            print(f"Affiliations: {', '.join(article.get('affiliations', [])) if article.get('affiliations') else 'Not available'}")
            print(f"URL: {article.get('url', 'Not available')}")
            if article.get("oa_pdf_path"):
                print(f"OA PDF Path: {article['oa_pdf_path']}")
            print("-" * 50)


if __name__ == "__main__":
    main()

