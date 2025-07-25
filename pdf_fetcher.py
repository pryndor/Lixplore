# pdf fetcher

import requests
import os

def fetch_pmc_pdf(pmcid_or_pmid, output_dir="downloads"):
    if pmcid_or_pmid.upper().startswith("PMC"):
        pmcid = pmcid_or_pmid
    else:
        # Convert PMID to PMCID using Entrez
        from Bio import Entrez
        Entrez.email = "your-email@example.com"
        Entrez.api_key = "your-api-key"

        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmcid_or_pmid)
        record = Entrez.read(handle)
        handle.close()

        try:
            pmcid = next(link["Id"] for link in record[0]["LinkSetDb"][0]["Link"])
            pmcid = "PMC" + pmcid
        except (IndexError, KeyError, StopIteration):
            print("❌ Could not find PMCID for the given PMID.")
            return

    url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf"

    try:
        os.makedirs(output_dir, exist_ok=True)
        file_path = os.path.join(output_dir, f"{pmcid}.pdf")
        r = requests.get(url)
        if r.ok and r.headers.get("Content-Type", "").startswith("application/pdf"):
            with open(file_path, "wb") as f:
                f.write(r.content)
            print(f"✅ PDF saved: {file_path}")
        else:
            print("❌ PDF not available or failed to download.")
    except Exception as e:
        print(f"⚠️ Error downloading PDF: {e}")

