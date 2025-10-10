# Testing file for scheduler

from apscheduler.schedulers.background import BackgroundScheduler
from datetime import datetime, timedelta
import time

# Import your PubMed search function
from api_clients.pubmed.pubmed_client import search_pubmed


# Define your scheduled job
def scheduled_pubmed_search():
    print(f"🔍 Running PubMed search at {datetime.now()}")
    try:
        results = search_pubmed("pharmacovigilance")  # Example query
        print("✅ PubMed search completed successfully!")
    except Exception as e:
        print(f"❌ Error during PubMed search: {e}")


# Initialize the scheduler
scheduler = BackgroundScheduler()

# Set your target run date and time
run_time = datetime(2025, 10, 10, 18, 17, 0)  # (year, month, day, hour, minute, second)

# Add job to run at that specific time
scheduler.add_job(scheduled_pubmed_search, 'date', run_date=run_time)

# Start the scheduler
scheduler.start()

print(f"🕒 Job scheduled for {run_time.strftime('%Y-%m-%d %H:%M:%S')}")
print("Scheduler running... Press Ctrl+C to exit.")

# Keep the program alive until manually stopped
try:
    while True:
        time.sleep(2)
except (KeyboardInterrupt, SystemExit):
    scheduler.shutdown()
    print("Scheduler stopped.")

