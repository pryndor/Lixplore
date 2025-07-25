import requests

def download_pdf(url, filename):
    try:
        response = requests.get(url)
        if response.status_code == 200 and 'application/pdf' in response.headers.get('Content-Type', ''):
            with open(filename, 'wb') as f:
                f.write(response.content)
            print(f"✅ Downloaded: {filename}")
        else:
            print("⚠️ Not a PDF or failed to fetch:", response.status_code)
    except Exception as e:
        print("❌ Error occurred:", e)

# Example Open Access PDF (You can try any direct PDF URL)
sample_url = "https://arxiv.org/pdf/2301.11982.pdf"
download_pdf(sample_url, "sample_article.pdf")

