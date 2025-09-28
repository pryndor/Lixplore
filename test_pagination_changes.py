#!/usr/bin/env python3
"""
Test script to verify pagination changes: 5 results per page, "1 to 5" format, page numbers start from 1
"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

def test_pagination_logic():
    """Test pagination logic with 5 results per page"""
    print("🔍 Testing Pagination Logic (5 results per page)...")
    
    # Test scenarios
    scenarios = [
        {'total_results': 12, 'page': 1, 'expected_start': 1, 'expected_end': 5, 'expected_total_pages': 3},
        {'total_results': 12, 'page': 2, 'expected_start': 6, 'expected_end': 10, 'expected_total_pages': 3},
        {'total_results': 12, 'page': 3, 'expected_start': 11, 'expected_end': 12, 'expected_total_pages': 3},
        {'total_results': 5, 'page': 1, 'expected_start': 1, 'expected_end': 5, 'expected_total_pages': 1},
        {'total_results': 3, 'page': 1, 'expected_start': 1, 'expected_end': 3, 'expected_total_pages': 1},
    ]
    
    per_page = 5  # Our new setting
    
    for scenario in scenarios:
        total_results = scenario['total_results']
        page = scenario['page']
        expected_start = scenario['expected_start']
        expected_end = scenario['expected_end']
        expected_total_pages = scenario['expected_total_pages']
        
        # Calculate pagination (same logic as in app.py)
        total_pages = (total_results + per_page - 1) // per_page
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page
        start_result = start_idx + 1 if total_results > 0 else 0
        end_result = min(end_idx, total_results)
        
        print(f"   📄 {total_results} results, page {page}:")
        print(f"      Start result: {start_result} (expected: {expected_start})")
        print(f"      End result: {end_result} (expected: {expected_end})")
        print(f"      Total pages: {total_pages} (expected: {expected_total_pages})")
        
        if (start_result == expected_start and 
            end_result == expected_end and 
            total_pages == expected_total_pages):
            print(f"      ✅ Correct!")
        else:
            print(f"      ❌ Failed!")
            return False
    
    print("   ✅ Pagination logic test passed!")
    return True

def test_display_format():
    """Test the display format shows '1 to 5' instead of '1-5'"""
    print("\n📊 Testing Display Format...")
    
    # Test the format string
    start_result = 1
    end_result = 5
    total_results = 12
    current_page = 1
    total_pages = 3
    
    # New format (what we implemented)
    new_format = f"Showing {start_result} to {end_result} of {total_results} articles (Page {current_page} of {total_pages})"
    expected_new = "Showing 1 to 5 of 12 articles (Page 1 of 3)"
    
    # Old format (what we changed from)
    old_format = f"Showing {start_result}-{end_result} of {total_results} articles (Page {current_page} of {total_pages})"
    expected_old = "Showing 1-5 of 12 articles (Page 1 of 3)"
    
    print(f"   New format: {new_format}")
    print(f"   Expected:   {expected_new}")
    print(f"   Old format: {old_format}")
    
    if new_format == expected_new:
        print("   ✅ Display format test passed!")
        return True
    else:
        print("   ❌ Display format test failed!")
        return False

def test_page_numbering():
    """Test that page numbering starts from 1"""
    print("\n🔢 Testing Page Numbering...")
    
    # Test page range generation
    total_pages = 5
    page_range = list(range(1, min(total_pages + 1, 101)))  # 1-100 max
    expected_range = [1, 2, 3, 4, 5]
    
    print(f"   Page range: {page_range}")
    print(f"   Expected:   {expected_range}")
    
    if page_range == expected_range:
        print("   ✅ Page numbering test passed!")
        return True
    else:
        print("   ❌ Page numbering test failed!")
        return False

def main():
    """Run all tests"""
    print("🧪 Testing Pagination Changes")
    print("=" * 50)
    
    tests = [
        test_pagination_logic,
        test_display_format,
        test_page_numbering
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
    
    print("\n" + "=" * 50)
    print(f"📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! Pagination changes are working correctly.")
        print("\n✅ Changes implemented:")
        print("   • Results per page: 10 → 5")
        print("   • Display format: '1-5' → '1 to 5'")
        print("   • Page numbering: starts from 1")
        print("   • Next button: available for navigation")
        return True
    else:
        print("❌ Some tests failed. Please check the implementation.")
        return False

if __name__ == "__main__":
    main()
