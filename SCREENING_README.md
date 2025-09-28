# 📋 Literature Search and Screening Tool

## Overview

The Literature Search and Screening Tool is a comprehensive addition to Lixplore that enables systematic literature reviews and meta-analyses. This tool provides a complete workflow for managing screening projects, coordinating multiple reviewers, and tracking progress through the screening process.

## 🌟 Key Features

### 1. **Screening Project Management**
- Create and manage multiple screening projects
- Define inclusion and exclusion criteria
- Configure dual screening requirements
- Set conflict resolution methods

### 2. **Multi-Reviewer Support**
- Add multiple reviewers to projects
- Assign different roles (Reviewer, Lead Reviewer, Admin)
- Track individual reviewer progress
- Support for dual screening workflows

### 3. **Article Screening Interface**
- Individual article screening with detailed view
- Batch screening for efficient processing
- Keyboard shortcuts for rapid decisions
- Include/Exclude/Maybe decision options with reasons

### 4. **Analytics and Reporting**
- Real-time screening progress tracking
- Inter-rater agreement calculations
- Reviewer performance metrics
- Visual progress charts and statistics

### 5. **Integration with Search Results**
- Direct import from Lixplore search results
- Create screening projects from existing searches
- Seamless workflow from search to screening

### 6. **Export Capabilities**
- Export screening results in CSV format
- Filter exports by screening status
- Include metadata and decision reasons

## 🚀 Getting Started

### 1. Access the Screening Dashboard
- Navigate to the main Lixplore interface
- Click on "📋 Screening" in the navigation menu
- Or visit `/screening` directly

### 2. Create Your First Screening Project
1. Click "➕ Create New Screening Project"
2. Fill in project details:
   - **Project Name**: Descriptive name for your review
   - **Description**: Brief overview of the project
   - **Lead Researcher**: Your name or primary contact
   - **Inclusion Criteria**: List what studies should be included
   - **Exclusion Criteria**: List what studies should be excluded
   - **Dual Screening**: Enable if you want two reviewers per article
   - **Conflict Resolution**: Choose how to handle disagreements

### 3. Add Reviewers
1. In your project dashboard, click "Add Reviewer"
2. Enter reviewer details:
   - Name and email
   - Role (Reviewer, Lead Reviewer, Admin)

### 4. Import Articles
You can add articles in several ways:

#### From Search Results:
1. Perform a search in Lixplore
2. Click "➕ Create Screening Project" in the results
3. Fill in project details and articles will be imported automatically

#### From Existing Searches:
1. In your project, click "➕ Add Articles"
2. Select "Import from Search Results"
3. Choose from your search history

#### Manual Entry:
1. Click "➕ Add Articles"
2. Select "Add Article Manually"
3. Enter article details manually

### 5. Start Screening
1. Go to your project dashboard
2. Click "Screen Article" on any pending article
3. Review the article details against your criteria
4. Make Include/Exclude/Maybe decisions
5. Add reasons and notes as needed

## 📊 Screening Workflow

### Individual Article Screening
1. **Review Article Details**: Title, abstract, authors, metadata
2. **Check Against Criteria**: Reference your inclusion/exclusion criteria
3. **Make Decision**: Include, Exclude, or Maybe
4. **Provide Reason**: Select from predefined reasons or add custom notes
5. **Continue**: Move to next article or return to project

### Batch Screening
1. Access batch screening from project dashboard
2. Review multiple articles in a streamlined interface
3. Use keyboard shortcuts (I=Include, E=Exclude, M=Maybe)
4. Make quick decisions with visual feedback
5. Save all decisions at once

### Conflict Resolution
When dual screening is enabled:
1. Articles with conflicting decisions are flagged
2. View conflicts in the project dashboard
3. Resolve through discussion, third reviewer, or consensus
4. Update final decisions

## 📈 Analytics and Progress Tracking

### Key Metrics
- **Total Articles**: Number of articles in the project
- **Screening Progress**: Percentage of articles screened
- **Decision Breakdown**: Count of included/excluded/pending articles
- **Inter-rater Agreement**: Percentage agreement between reviewers
- **Completion Rate**: Overall project progress

### Reviewer Performance
- Individual decision counts
- Inclusion rates by reviewer
- Screening velocity (decisions per day)
- Progress comparison between reviewers

### Visual Reports
- Progress bar showing screening status
- Daily activity charts
- Reviewer performance tables
- Export-ready summary statistics

## 📥 Export and Reporting

### Available Export Formats
- **CSV**: Comma-separated values for spreadsheet analysis
- **Filtered Exports**: Export by status (included, excluded, conflicts)
- **Complete Dataset**: All articles with metadata and decisions

### Export Options
1. **All Results**: Complete dataset with all articles
2. **Included Articles**: Only articles marked for inclusion
3. **Excluded Articles**: Articles excluded with reasons
4. **Conflicts**: Articles requiring resolution

## 🔧 Advanced Features

### Dual Screening
- Enable for systematic reviews requiring two independent reviewers
- Automatic conflict detection when reviewers disagree
- Configurable conflict resolution workflows

### Keyboard Shortcuts
In batch screening mode:
- `I` - Mark as Include
- `E` - Mark as Exclude  
- `M` - Mark as Maybe

### Project Status Management
- **Active**: Project is ongoing
- **Completed**: Screening is finished
- **Paused**: Temporarily suspended

## 🛠️ Technical Details

### Database Schema
The screening functionality adds several new tables:
- `screening_projects`: Project metadata and settings
- `screening_reviewers`: Reviewer information and roles
- `screening_articles`: Articles linked to projects
- `screening_decisions`: Individual screening decisions
- `screening_conflicts`: Conflict tracking and resolution

### Integration Points
- Seamless integration with existing Lixplore search functionality
- Compatible with all supported databases (PubMed, CrossRef, etc.)
- Maintains existing user interface patterns and styling

## 🎯 Best Practices

### Project Setup
1. **Clear Criteria**: Define specific, measurable inclusion/exclusion criteria
2. **Pilot Testing**: Screen a small sample to refine criteria
3. **Reviewer Training**: Ensure all reviewers understand the criteria
4. **Regular Meetings**: Schedule progress reviews and conflict resolution

### Screening Process
1. **Consistent Application**: Apply criteria uniformly across all articles
2. **Document Decisions**: Use reason codes and notes for transparency
3. **Regular Backups**: Export progress regularly
4. **Quality Checks**: Monitor inter-rater agreement

### Team Coordination
1. **Role Assignment**: Clearly define reviewer roles and responsibilities
2. **Communication**: Use project notes for team communication
3. **Progress Monitoring**: Regular check-ins on screening progress
4. **Conflict Resolution**: Establish clear procedures for disagreements

## 🆘 Troubleshooting

### Common Issues
1. **Articles Not Importing**: Check search history and try manual import
2. **Reviewer Access**: Verify reviewer email addresses and roles
3. **Slow Performance**: Consider batch screening for large datasets
4. **Export Issues**: Ensure sufficient articles are screened before export

### Support
For technical support or feature requests:
1. Check the main Lixplore documentation
2. Contact the development team
3. Submit issues through the project repository

## 🔮 Future Enhancements

Planned features for future releases:
- **Bulk Import**: CSV/BibTeX file import
- **Advanced Analytics**: Machine learning screening assistance
- **API Integration**: Connect with reference management tools
- **Mobile Interface**: Responsive design for mobile screening
- **Collaboration Tools**: Real-time collaboration features

---

**Ready to start your systematic literature review?** 
Navigate to the Screening Dashboard and create your first project!
