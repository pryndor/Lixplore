import requests

pmid = "31801969"  # example PMID
url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
params = {"ids": pmid, "format": "json"}
response = requests.get(url, params=params)
data = response.json()
print(data)

