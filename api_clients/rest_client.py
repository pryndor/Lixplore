# This file handling Rest API for clients

import requests

def fetch_paginated_data(base_url, params=None, headers=None, max_pages=None):
    """Generic REST API pagination handler"""
    all_results = []
    page = 1

    while True:
        params = params or {}
        params["page"] = page
        response = requests.get(base_url, params=params, headers=headers)
        if response.status_code != 200:
            break

        data = response.json()
        results = data.get("results", data)
        if not results:
            break

        all_results.extend(results)

        # GitHub-like pagination check
        if "next" not in response.links:
            break

        page += 1
        if max_pages and page > max_pages:
            break

    return all_results



