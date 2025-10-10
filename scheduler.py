# This is a scheduler script

from apscheduler.schedulers.background import BackgroundScheduler
from datetime import datetime
import time
from api_clients.pubmed.pubmed_client import search_pubmed

def automated_search():
    keywords = ["pharmacovigilance", "drug safety", "signal detection"]
    print(f"[{datetime.now()}] Running automated PubMed search...")
    all_results = []

    for keyword in keywords:
        print(f"🔍 Searching PubMed for: {keyword}")
        results = search_pubmed(keyword)
        if results:
            all_results.extend(results)
        print(f"✅ Retrieved {len(results) if results else 0} results for '{keyword}'")

    print(f"🎯 Total combined results: {len(all_results)}")
    # store_results_in_db(all_results)  # Optional if you have DB function

scheduler = BackgroundScheduler()
scheduler.add_job(automated_search, 'interval', hours=6)  # Run every 6 hours
scheduler.start()

print("🕒 Scheduler started. Waiting for jobs...")

try:
    while True:
        time.sleep(2)
except (KeyboardInterrupt, SystemExit):
    scheduler.shutdown()
    print("🛑 Scheduler stopped.")

