# 🔍 Enhanced Search History Data Management

## Overview

The Search History feature now includes comprehensive data management with Local Search and Global Search capabilities, providing users with powerful filtering and search options within their search history.

## ✅ **Enhanced Features**

### 🏠 **Local Search**
- **Scope**: Searches only in the query text field
- **Use Case**: Find searches by specific query terms
- **Example**: Searching for "COVID" will find all searches containing "COVID" in the query
- **Visual Indicator**: 🏠 Query Match badge

### 🌐 **Global Search**
- **Scope**: Searches across query, source, and timestamp fields
- **Use Case**: Find searches by any criteria (query, database source, or date)
- **Example**: Searching for "pubmed" will find all PubMed searches
- **Visual Indicator**: 🌐 Global Match badge with field details

### 📊 **Enhanced Data Structure**

#### **Database Enhancements**
```python
# Enhanced get_archive() function
def get_archive():
    """Retrieve search history with article counts and enhanced data"""
    # Returns dictionary format with:
    # - id, source, query, filters, timestamp
    # - article_count (number of results)
    # - search_type, display_type
    # - is_local_match, is_global_match
    # - matched_fields (for highlighting)
```

#### **New Database Functions**
```python
def get_search_history_by_type(search_type, search_query, source_filter):
    """Categorized search with filtering by type"""

def get_matched_fields(row, search_query, search_type):
    """Determine which fields matched for highlighting"""
```

### 🎨 **Visual Enhancements**

#### **Search Term Highlighting**
- **Local Search**: Highlights matching terms in query field
- **Global Search**: Highlights matching terms in all fields (query, source, timestamp)
- **Implementation**: Uses `<mark>` tags with yellow background
- **Case Insensitive**: Matches regardless of case

#### **Match Indicators**
- **Local Match Badge**: 🏠 Green badge for query-only matches
- **Global Match Badge**: 🌐 Blue badge for multi-field matches
- **Field Tags**: Small tags showing which fields matched (Query, Source, Date)

#### **Enhanced Cards**
- **Article Count**: Shows number of results for each search
- **Search Result Styling**: Special border and gradient for search results
- **Improved Actions**: "View X Results" button shows exact count

## 🔧 **Technical Implementation**

### **Backend Changes (app.py)**

#### **Enhanced Search History Route**
```python
@app.route("/search-history")
def search_history():
    # Get filtered history using new database functions
    filtered_history = get_search_history_by_type(
        search_type=search_type,
        search_query=search_query,
        source_filter=source_filter
    )
    
    # Add match information for highlighting
    for item in filtered_history:
        item['display_type'] = search_type
        item['is_local_match'] = search_query and search_query.lower() in item.get('query', '').lower()
        item['is_global_match'] = search_query and (
            search_query.lower() in item.get('query', '').lower() or
            search_query.lower() in item.get('source', '').lower() or
            search_query.lower() in str(item.get('timestamp', '')).lower()
        )
```

### **Database Changes (database.py)**

#### **Enhanced Archive Function**
- **Article Count**: JOIN with articles table to count results
- **Dictionary Format**: Returns structured data instead of tuples
- **Search Type**: Adds metadata for search categorization

#### **New Filtering Function**
- **Type-Based Search**: Separate logic for local vs global
- **Field Matching**: Tracks which fields matched the search
- **Source Filtering**: Combined with search type filtering

### **Frontend Changes (search_history.html)**

#### **Enhanced Search Controls**
```html
<div class="search-stats">
    <p><strong>🔍 Searching for:</strong> "{{ search_query }}"</p>
    <p><strong>📊 Search Mode:</strong> {{ search_type.title() }} Search</p>
</div>
```

#### **Dynamic Highlighting**
```html
{% if search_query and search.is_local_match %}
    {{ search.query | replace(search_query, '<mark>' + search_query + '</mark>') | safe }}
{% else %}
    {{ search.query }}
{% endif %}
```

#### **Match Information Display**
```html
<div class="match-info">
    {% if search_type == 'local' %}
        <span class="match-badge local-match">🏠 Query Match</span>
    {% else %}
        <span class="match-badge global-match">🌐 Global Match</span>
        <div class="match-details">
            <span class="field-match">Query</span>
            <span class="field-match">Source</span>
            <span class="field-match">Date</span>
        </div>
    {% endif %}
</div>
```

## 🎯 **User Experience**

### **Search Workflow**
1. **Access**: Navigate to "🔍 Search History"
2. **Choose Mode**: Select Local or Global search
3. **Enter Query**: Type search term
4. **Filter Source**: Optionally filter by database
5. **View Results**: See highlighted matches with indicators

### **Visual Feedback**
- **Search Type Info**: Clear explanation of current search mode
- **Search Stats**: Shows current search term and mode
- **Highlighted Terms**: Yellow highlighting on matched text
- **Match Badges**: Color-coded indicators for match type
- **Field Tags**: Shows which fields matched in global search
- **Article Counts**: Displays result count for each search

### **Data Organization**
- **Alphabetical Sorting**: All results sorted by query text
- **Pagination**: 10 results per page with enhanced controls
- **Source Filtering**: Filter by specific databases
- **Article Counts**: Shows productivity of each search

## 📊 **Data Categories**

### **Local Search Results**
- **Criteria**: Query field contains search term
- **Display**: 🏠 Query Match badge
- **Highlighting**: Only query field highlighted
- **Use Case**: Finding specific research topics

### **Global Search Results**
- **Criteria**: Any field contains search term
- **Display**: 🌐 Global Match badge + field tags
- **Highlighting**: All matching fields highlighted
- **Use Case**: Finding searches by database, date, or any criteria

### **Combined Results**
- **Article Counts**: Shows research productivity
- **Timestamps**: Chronological context
- **Source Information**: Database diversity
- **Filter Options**: Refined searching

## 🎨 **Styling Features**

### **Color Scheme**
- **Local Match**: Green (#d4edda) - focused search
- **Global Match**: Blue (#cce5ff) - comprehensive search
- **Highlighting**: Yellow (#fff3cd) - search terms
- **Field Tags**: Gray (#e9ecef) - field indicators

### **Interactive Elements**
- **Hover Effects**: Cards lift on hover
- **Search Result Cards**: Special border and gradient
- **Responsive Design**: Works on all screen sizes
- **Touch Friendly**: Adequate button sizes

## 🔮 **Advanced Features**

### **Smart Highlighting**
- **Case Insensitive**: Matches regardless of case
- **Partial Matches**: Highlights partial word matches
- **Multiple Fields**: Highlights in query, source, and timestamp
- **Safe Rendering**: Uses Jinja2 safe filter for HTML

### **Match Detection**
- **Field-Specific**: Knows which fields matched
- **Visual Indicators**: Shows match scope clearly
- **Search Context**: Provides context for matches
- **Performance Optimized**: Efficient string matching

## 📋 **Summary**

### ✅ **Local Search Features**
- **🏠 Query-Only Search**: Searches only in query text field
- **🎯 Focused Results**: Finds specific research topics
- **💡 Query Highlighting**: Highlights matching terms in queries
- **📊 Article Counts**: Shows result count for each search

### ✅ **Global Search Features**
- **🌐 Multi-Field Search**: Searches across query, source, and timestamp
- **🔍 Comprehensive Results**: Finds searches by any criteria
- **🎨 Full Highlighting**: Highlights matching terms in all fields
- **🏷️ Field Tags**: Shows which fields matched (Query, Source, Date)

### ✅ **Enhanced Data Management**
- **📊 Article Counts**: Shows productivity of each search
- **📝 Alphabetical Sorting**: Organized result display
- **🔧 Source Filtering**: Filter by specific databases
- **📄 Enhanced Pagination**: 10 results per page with circular navigation

### ✅ **Visual Improvements**
- **🎨 Search Highlighting**: Yellow highlighting for matched terms
- **🏷️ Match Badges**: Color-coded indicators for match type
- **📱 Responsive Design**: Works on all devices
- **✨ Interactive Effects**: Hover animations and visual feedback

The enhanced Search History now provides a comprehensive data management system that allows users to efficiently search, filter, and organize their research history with clear visual indicators and powerful search capabilities!
