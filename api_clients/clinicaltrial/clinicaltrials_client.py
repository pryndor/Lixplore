# clinicaltrials_client.py

# clinicaltrials_client.py
import requests

BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

def search_clinical_trials(keyword, status=None, page_size=5, max_pages=2):
    """
    Search ClinicalTrials.gov for studies.
    
    :param keyword: Search term (e.g., "paracetamol")
    :param status: Recruitment status filter (e.g., "RECRUITING")
    :param page_size: Number of results per page
    :param max_pages: Max number of pages to fetch
    :return: List of study details
    """
    params = {
        "query.cond": keyword,
        "pageSize": page_size
    }
    
    if status:
        params["filter.overallStatus"] = status

    studies_list = []
    page_count = 0
    next_page_token = None

    while page_count < max_pages:
        if next_page_token:
            params["pageToken"] = next_page_token
        
        response = requests.get(BASE_URL, params=params)
        
        if response.status_code != 200:
            print("Error:", response.status_code, response.text)
            break
        
        data = response.json()
        studies = data.get("studies", [])
        
        for study in studies:
            protocol = study.get("protocolSection", {})
            identification = protocol.get("identificationModule", {})
            status_module = protocol.get("statusModule", {})
            conditions_module = protocol.get("conditionsModule", {})
            interventions_module = protocol.get("armsInterventionsModule", {})
            description_module = protocol.get("descriptionModule", {})
            design_module = protocol.get("designModule", {})

            nct_id = identification.get("nctId", "N/A")
            title = identification.get("briefTitle", "N/A")
            overall_status = status_module.get("overallStatus", "N/A")
            conditions = conditions_module.get("conditions", [])
            interventions = interventions_module.get("interventions", [])
            study_type = design_module.get("studyType", "N/A")
            phase = design_module.get("phase", "N/A")
            start_date = status_module.get("startDateStruct", {}).get("date", "N/A")
            completion_date = status_module.get("completionDateStruct", {}).get("date", "N/A")
            summary = description_module.get("briefSummary", "N/A")
            
            url = f"https://clinicaltrials.gov/study/{nct_id}" if nct_id != "N/A" else "N/A"

            studies_list.append({
                "nct_id": nct_id,
                "title": title,
                "status": overall_status,
                "conditions": conditions,
                "interventions": interventions,
                "study_type": study_type,
                "phase": phase,
                "start_date": start_date,
                "completion_date": completion_date,
                "summary": summary,
                "url": url
            })
        
        next_page_token = data.get("nextPageToken")
        if not next_page_token:
            break
        
        page_count += 1

    return studies_list

