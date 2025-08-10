ATTACH 'lixplore.db' AS extra;

INSERT INTO articles (
    search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path, created_at
)
SELECT 
    search_id, title, authors, doi, pmid, url, abstract, affiliations, oa_pdf_path, CURRENT_TIMESTAMP
FROM extra.articles;

-- Repeat for other tables...

DETACH extra;

