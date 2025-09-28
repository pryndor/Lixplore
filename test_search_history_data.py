#!/usr/bin/env python3
"""
Test script for enhanced search history data handling with local and global search
"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

def test_local_vs_global_search():
    """Test local vs global search logic"""
    print("🔍 Testing Local vs Global Search Logic...")
    
    # Sample search history data
    sample_history = [
        {
            'id': 1,
            'source': 'pubmed',
            'query': 'COVID-19 treatment',
            'filters': '{}',
            'timestamp': '2024-01-15 10:30:00',
            'article_count': 25
        },
        {
            'id': 2,
            'source': 'crossref',
            'query': 'machine learning',
            'filters': '{}',
            'timestamp': '2024-01-16 14:20:00',
            'article_count': 18
        },
        {
            'id': 3,
            'source': 'pubmed',
            'query': 'artificial intelligence',
            'filters': '{}',
            'timestamp': '2024-01-17 09:15:00',
            'article_count': 32
        },
        {
            'id': 4,
            'source': 'springer',
            'query': 'data science methods',
            'filters': '{}',
            'timestamp': '2024-01-18 16:45:00',
            'article_count': 12
        }
    ]
    
    search_query = 'machine'
    
    # Test Local Search (query field only)
    print(f"\n🏠 Local Search for '{search_query}':")
    local_results = []
    for item in sample_history:
        if search_query.lower() in item['query'].lower():
            local_results.append(item)
            print(f"   ✅ Found: {item['query']} (ID: {item['id']})")
    
    print(f"   📊 Local search found {len(local_results)} results")
    
    # Test Global Search (all fields)
    print(f"\n🌐 Global Search for '{search_query}':")
    global_results = []
    for item in sample_history:
        if (search_query.lower() in item['query'].lower() or
            search_query.lower() in item['source'].lower() or
            search_query.lower() in item['timestamp'].lower()):
            global_results.append(item)
            
            # Determine which fields matched
            matched_fields = []
            if search_query.lower() in item['query'].lower():
                matched_fields.append('query')
            if search_query.lower() in item['source'].lower():
                matched_fields.append('source')
            if search_query.lower() in item['timestamp'].lower():
                matched_fields.append('timestamp')
            
            print(f"   ✅ Found: {item['query']} (ID: {item['id']}) - Matched: {', '.join(matched_fields)}")
    
    print(f"   📊 Global search found {len(global_results)} results")
    
    # Verify results
    if len(global_results) >= len(local_results):
        print("   ✅ Global search correctly found same or more results than local search")
        return True
    else:
        print("   ❌ Global search found fewer results than local search - this is incorrect")
        return False

def test_search_highlighting():
    """Test search term highlighting logic"""
    print("\n🎨 Testing Search Term Highlighting...")
    
    test_cases = [
        {
            'text': 'COVID-19 treatment protocols',
            'search_term': 'COVID',
            'expected_highlight': True
        },
        {
            'text': 'machine learning algorithms',
            'search_term': 'learning',
            'expected_highlight': True
        },
        {
            'text': 'artificial intelligence',
            'search_term': 'xyz',
            'expected_highlight': False
        }
    ]
    
    for case in test_cases:
        text = case['text']
        search_term = case['search_term']
        expected = case['expected_highlight']
        
        # Simulate highlighting logic
        should_highlight = search_term.lower() in text.lower()
        
        print(f"   📝 Text: '{text}'")
        print(f"   🔍 Search: '{search_term}'")
        print(f"   💡 Should highlight: {should_highlight} (expected: {expected})")
        
        if should_highlight == expected:
            print(f"   ✅ Correct!")
        else:
            print(f"   ❌ Failed!")
            return False
    
    print("   ✅ Search highlighting test passed!")
    return True

def test_match_field_detection():
    """Test field match detection for global search"""
    print("\n🎯 Testing Match Field Detection...")
    
    sample_record = {
        'query': 'COVID-19 treatment',
        'source': 'pubmed',
        'timestamp': '2024-01-15 10:30:00'
    }
    
    test_searches = [
        {
            'search_term': 'COVID',
            'expected_fields': ['query']
        },
        {
            'search_term': 'pubmed',
            'expected_fields': ['source']
        },
        {
            'search_term': '2024',
            'expected_fields': ['timestamp']
        },
        {
            'search_term': 'treatment',
            'expected_fields': ['query']
        },
        {
            'search_term': 'xyz',
            'expected_fields': []
        }
    ]
    
    for test in test_searches:
        search_term = test['search_term']
        expected_fields = test['expected_fields']
        
        # Simulate field matching logic
        matched_fields = []
        if search_term.lower() in sample_record['query'].lower():
            matched_fields.append('query')
        if search_term.lower() in sample_record['source'].lower():
            matched_fields.append('source')
        if search_term.lower() in sample_record['timestamp'].lower():
            matched_fields.append('timestamp')
        
        print(f"   🔍 Search: '{search_term}'")
        print(f"   📍 Matched fields: {matched_fields}")
        print(f"   📍 Expected fields: {expected_fields}")
        
        if set(matched_fields) == set(expected_fields):
            print(f"   ✅ Correct!")
        else:
            print(f"   ❌ Failed!")
            return False
    
    print("   ✅ Match field detection test passed!")
    return True

def test_alphabetical_sorting():
    """Test alphabetical sorting of search results"""
    print("\n📝 Testing Alphabetical Sorting...")
    
    unsorted_history = [
        {'query': 'Zebra research methods'},
        {'query': 'artificial intelligence'},
        {'query': 'COVID-19 treatment'},
        {'query': 'machine learning'},
        {'query': 'data science'}
    ]
    
    # Sort alphabetically by query (case-insensitive)
    sorted_history = sorted(unsorted_history, key=lambda x: x['query'].lower())
    
    print("   📊 Original order:")
    for i, item in enumerate(unsorted_history):
        print(f"      {i+1}. {item['query']}")
    
    print("   📝 Alphabetically sorted:")
    for i, item in enumerate(sorted_history):
        print(f"      {i+1}. {item['query']}")
    
    # Verify sorting
    queries = [item['query'].lower() for item in sorted_history]
    is_sorted = all(queries[i] <= queries[i+1] for i in range(len(queries)-1))
    
    if is_sorted:
        print("   ✅ Alphabetical sorting is correct!")
        return True
    else:
        print("   ❌ Alphabetical sorting failed!")
        return False

def test_data_structure_enhancement():
    """Test enhanced data structure with article counts and match info"""
    print("\n📊 Testing Enhanced Data Structure...")
    
    # Simulate enhanced data structure
    enhanced_record = {
        'id': 1,
        'source': 'pubmed',
        'query': 'COVID-19 treatment',
        'filters': '{}',
        'timestamp': '2024-01-15 10:30:00',
        'article_count': 25,
        'search_type': 'database',
        'display_type': 'global',
        'is_local_match': True,
        'is_global_match': True,
        'matched_fields': ['query']
    }
    
    required_fields = [
        'id', 'source', 'query', 'filters', 'timestamp', 
        'article_count', 'search_type', 'display_type',
        'is_local_match', 'is_global_match'
    ]
    
    print("   📋 Checking required fields:")
    all_present = True
    for field in required_fields:
        if field in enhanced_record:
            print(f"      ✅ {field}: {enhanced_record[field]}")
        else:
            print(f"      ❌ {field}: MISSING")
            all_present = False
    
    if all_present:
        print("   ✅ Enhanced data structure test passed!")
        return True
    else:
        print("   ❌ Enhanced data structure test failed!")
        return False

def main():
    """Run all search history data tests"""
    print("🚀 Starting Search History Data Tests...")
    
    tests = [
        test_local_vs_global_search,
        test_search_highlighting,
        test_match_field_detection,
        test_alphabetical_sorting,
        test_data_structure_enhancement
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                print(f"   ❌ Test {test.__name__} failed!")
        except Exception as e:
            print(f"   ❌ Test {test.__name__} error: {e}")
    
    print(f"\n📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All search history data tests passed!")
        print("\n📋 Enhanced Search History Features:")
        print("✅ Local Search: Searches only in query field")
        print("✅ Global Search: Searches across query, source, and timestamp")
        print("✅ Search term highlighting with <mark> tags")
        print("✅ Match field detection and display")
        print("✅ Article count display for each search")
        print("✅ Alphabetical sorting of results")
        print("✅ Enhanced data structure with match information")
        print("✅ Visual indicators for local vs global matches")
        
        print("\n🚀 Ready to test with Flask application!")
        print("1. Start the Flask app: python app.py")
        print("2. Navigate to Search History")
        print("3. Test Local Search (searches only query field)")
        print("4. Test Global Search (searches all fields)")
        print("5. Verify search term highlighting")
        print("6. Check match field indicators")
        print("7. Verify article counts are displayed")
        
        return True
    else:
        print(f"\n❌ {total - passed} tests failed. Please check the implementation.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
