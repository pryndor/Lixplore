#!/usr/bin/env python3
"""
Test script for the Literature Screening functionality
"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

from database.database import (
    init_db, create_screening_project, add_reviewer_to_project,
    add_articles_to_screening_project, get_screening_projects,
    get_screening_statistics
)

def test_screening_functionality():
    """Test the basic screening functionality"""
    print("🧪 Testing Literature Screening Functionality...")
    
    # Initialize database
    print("1. Initializing database...")
    init_db()
    print("   ✅ Database initialized")
    
    # Create a test screening project
    print("2. Creating test screening project...")
    project_id = create_screening_project(
        name="Test COVID-19 Systematic Review",
        description="A test screening project for COVID-19 treatment studies",
        inclusion_criteria="- RCT studies\n- COVID-19 patients\n- Published 2020-2024",
        exclusion_criteria="- Case reports\n- Animal studies\n- Non-English",
        created_by="Test Researcher",
        require_dual_screening=True,
        conflict_resolution_method="discussion"
    )
    
    if project_id:
        print(f"   ✅ Project created with ID: {project_id}")
    else:
        print("   ❌ Failed to create project")
        return False
    
    # Add reviewers
    print("3. Adding reviewers...")
    reviewer1_id = add_reviewer_to_project(project_id, "Dr. Alice Smith", "alice@example.com", "reviewer")
    reviewer2_id = add_reviewer_to_project(project_id, "Dr. Bob Johnson", "bob@example.com", "reviewer")
    
    if reviewer1_id and reviewer2_id:
        print(f"   ✅ Reviewers added: {reviewer1_id}, {reviewer2_id}")
    else:
        print("   ❌ Failed to add reviewers")
        return False
    
    # Add test articles
    print("4. Adding test articles...")
    test_articles = [
        {
            'title': 'Efficacy of Remdesivir in COVID-19 Patients: A Randomized Controlled Trial',
            'authors': 'Smith A, Johnson B, Williams C',
            'abstract': 'This study evaluates the efficacy of remdesivir in treating COVID-19 patients...',
            'doi': '10.1000/test1',
            'pmid': '12345678',
            'url': 'https://example.com/article1',
            'publication_year': 2021,
            'journal': 'New England Journal of Medicine'
        },
        {
            'title': 'Dexamethasone Treatment in Severe COVID-19: A Meta-Analysis',
            'authors': 'Brown D, Davis E, Miller F',
            'abstract': 'Meta-analysis of dexamethasone treatment outcomes in severe COVID-19 cases...',
            'doi': '10.1000/test2',
            'pmid': '87654321',
            'url': 'https://example.com/article2',
            'publication_year': 2022,
            'journal': 'The Lancet'
        },
        {
            'title': 'COVID-19 Vaccine Efficacy in Elderly Populations',
            'authors': 'Wilson G, Taylor H, Anderson I',
            'abstract': 'Study of COVID-19 vaccine efficacy in elderly populations...',
            'doi': '10.1000/test3',
            'pmid': '11223344',
            'url': 'https://example.com/article3',
            'publication_year': 2023,
            'journal': 'JAMA'
        }
    ]
    
    articles_added = add_articles_to_screening_project(project_id, test_articles)
    if articles_added > 0:
        print(f"   ✅ Added {articles_added} test articles")
    else:
        print("   ❌ Failed to add articles")
        return False
    
    # Get project statistics
    print("5. Getting project statistics...")
    statistics = get_screening_statistics(project_id)
    if statistics:
        print(f"   ✅ Statistics retrieved:")
        print(f"      - Total articles: {statistics['total_articles']}")
        print(f"      - Status counts: {statistics['status_counts']}")
        print(f"      - Reviewer progress: {len(statistics['reviewer_progress'])} reviewers")
    else:
        print("   ❌ Failed to get statistics")
        return False
    
    # List all projects
    print("6. Listing all screening projects...")
    projects = get_screening_projects()
    if projects:
        print(f"   ✅ Found {len(projects)} screening projects")
        for project in projects:
            print(f"      - {project['name']} (Status: {project['status']})")
    else:
        print("   ❌ No projects found")
        return False
    
    print("\n🎉 All screening functionality tests passed!")
    print("\n📋 Literature Screening Tool is ready to use!")
    print("\nNext steps:")
    print("1. Start the Flask application: python app.py")
    print("2. Navigate to http://localhost:5000")
    print("3. Click on '📋 Screening' in the navigation")
    print("4. Create your first screening project")
    print("5. Import articles from search results or add manually")
    print("6. Start screening articles with your team")
    
    return True

if __name__ == "__main__":
    success = test_screening_functionality()
    if not success:
        print("\n❌ Some tests failed. Please check the error messages above.")
        sys.exit(1)
    else:
        print("\n✅ All tests completed successfully!")
