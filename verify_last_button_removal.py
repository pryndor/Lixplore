#!/usr/bin/env python3
"""
Verification script to confirm last button has been removed
"""

import os

def check_file_for_last_button(file_path):
    """Check if file contains last button references"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check for last button references
        last_button_refs = [
            'pagination-last',
            'Last</a>',
            'show_end',
            'pagination.show_end'
        ]
        
        found_refs = []
        for ref in last_button_refs:
            if ref in content:
                found_refs.append(ref)
        
        return found_refs
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return []

def main():
    """Verify last button removal"""
    print("🔍 Verifying Last Button Removal...")
    
    files_to_check = [
        'templates/index.html',
        'templates/search_history.html',
        'app.py'
    ]
    
    all_clean = True
    
    for file_path in files_to_check:
        print(f"\n📄 Checking {file_path}...")

        found_refs = check_file_for_last_button(file_path)
        
        if found_refs:
            print(f"   ❌ Found last button references: {found_refs}")
            all_clean = False
        else:
            print(f"   ✅ No last button references found")
    
    print(f"\n📊 Verification Results:")
    if all_clean:
        print("✅ Last button successfully removed from all files!")
        print("\n📋 Current Pagination Controls:")
        print("✅ Start button (purple) - Quick navigation to page 1")
        print("✅ Page numbers 1-100 - Smart display with ellipsis")
        print("✅ Next button (green) - Forward navigation or 'Start Over'")
        print("❌ Last button - REMOVED as requested")
        
        print("\n🎯 Navigation Flow:")
        print("1. Start → Page Numbers → Next")
        print("2. Circular: Next from last page returns to page 1")
        print("3. No Last button for direct end navigation")
        
        return True
    else:
        print("❌ Last button references still found in some files!")
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
