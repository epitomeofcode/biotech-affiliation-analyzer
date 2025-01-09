"""
PubMed Paper Analyzer
Fetches and analyzes research papers with industry affiliations.
"""
import argparse
import sys
import logging
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime
from Bio import Entrez
import pandas as pd
import re
from pathlib import Path

@dataclass
class AuthorInfo:
    """Structured container for author information"""
    name: str
    affiliation: str
    email: Optional[str] = None
    is_corresponding: bool = False

@dataclass
class PaperData:
    """Structured container for paper information"""
    pubmed_id: str
    title: str
    pub_date: str
    industry_authors: List[AuthorInfo]
    
    def to_dict(self) -> Dict:
        """Convert paper data to dictionary format for CSV export"""
        return {
            'PubmedID': self.pubmed_id,
            'Title': self.title,
            'Publication_Date': self.pub_date,
            'Non_academic_Authors': '; '.join(
                [author.name for author in self.industry_authors]
            ),
            'Company_Affiliations': '; '.join(
                set(author.affiliation for author in self.industry_authors)
            ),
            'Corresponding_Author_Email': next(
                (author.email for author in self.industry_authors 
                 if author.is_corresponding), 
                ''
            )
        }

class IndustryAffiliationAnalyzer:
    """Analyzes affiliations to identify industry-affiliated authors"""
    
    def __init__(self):
        # Custom curated lists for affiliation analysis
        self._industry_indicators = {
            'core_terms': {'inc', 'corp', 'ltd', 'llc', 'company', 'industries'},
            'biotech_terms': {'therapeutics', 'biotechnology', 'bioscience', 
                            'pharmaceuticals', 'biopharma', 'laboratories'},
            'divisions': {'research', 'development', 'discovery', 'clinical', 
                        'regulatory'}
        }
        self._academic_indicators = {
            'institutions': {'university', 'college', 'institute', 'school'},
            'medical': {'hospital', 'clinic', 'medical center', 'health system'},
            'research': {'national laboratory', 'research center', 'institute of'}
        }
        
    def analyze_affiliation(self, affiliation: str) -> Tuple[bool, str]:
        """
        Analyzes an affiliation string to determine if it's an industry affiliation
        Returns: (is_industry, cleaned_affiliation_name)
        """
        if not affiliation:
            return False, ''
            
        text = affiliation.lower()
        
        # Quick reject for obvious academic institutions
        for academic_type in self._academic_indicators.values():
            if any(term in text for term in academic_type):
                return False, ''
        
        # Check for industry indicators
        industry_score = 0
        for category in self._industry_indicators.values():
            if any(term in text for term in category):
                industry_score += 1
                
        # Extract organization name using pattern matching
        org_pattern = r'([A-Z][A-Za-z0-9\s&]+(?:Inc\.|Corp\.|Ltd\.|LLC|GmbH)?)'
        org_match = re.search(org_pattern, affiliation)
        org_name = org_match.group(1).strip() if org_match else affiliation
        
        return industry_score >= 2, org_name

class PubMedResearchFetcher:
    """Handles PubMed API interaction and paper processing"""
    
    def __init__(self, contact_email: str):
        self.contact_email = contact_email
        Entrez.email = contact_email
        self._affiliation_analyzer = IndustryAffiliationAnalyzer()
        self._logger = logging.getLogger(__name__)
        
    def _extract_email(self, text: str) -> Optional[str]:
        """Extract email address from text using custom pattern matching"""
        email_pattern = r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}'
        match = re.search(email_pattern, text)
        return match.group(0) if match else None
        
    def _parse_publication_date(self, date_str: str) -> str:
        """Convert PubMed date format to consistent output format"""
        try:
            # Handle various date formats
            date_patterns = [
                '%Y %b %d',
                '%Y %b',
                '%Y'
            ]
            
            for pattern in date_patterns:
                try:
                    date_obj = datetime.strptime(date_str, pattern)
                    return date_obj.strftime('%Y-%m-%d')
                except ValueError:
                    continue
                    
            return date_str  # Return original if no pattern matches
            
        except Exception as e:
            self._logger.warning(f"Date parsing error: {e}")
            return date_str
            
    def _process_record(self, record: Dict) -> Optional[PaperData]:
        """Process a single PubMed record into structured paper data"""
        try:
            # Extract and process authors
            industry_authors = []
            for author in record.get('FAU', []):
                author_data = record.get('AD', [])
                if not author_data:
                    continue
                    
                # Process each affiliation
                for affiliation in author_data:
                    is_industry, org_name = (
                        self._affiliation_analyzer.analyze_affiliation(affiliation)
                    )
                    
                    if is_industry:
                        # Check if corresponding author
                        is_corresponding = (
                            'corresponding author' in affiliation.lower()
                        )
                        email = (
                            self._extract_email(affiliation) 
                            if is_corresponding 
                            else None
                        )
                        
                        industry_authors.append(AuthorInfo(
                            name=author,
                            affiliation=org_name,
                            email=email,
                            is_corresponding=is_corresponding
                        ))
                        
            if not industry_authors:
                return None
                
            return PaperData(
                pubmed_id=record.get('PMID', ''),
                title=record.get('TI', ''),
                pub_date=self._parse_publication_date(record.get('DP', '')),
                industry_authors=industry_authors
            )
            
        except Exception as e:
            self._logger.error(f"Record processing error: {e}")
            return None
            
    def fetch_papers(self, query: str) -> List[PaperData]:
        """
        Fetch and process papers from PubMed matching the query
        """
        self._logger.info(f"Initiating PubMed search: {query}")
        
        try:
            # Perform initial search
            search_handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=1000,
                usehistory="y"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if not search_results['Count']:
                self._logger.info("No matching papers found")
                return []
                
            # Fetch full records using history
            fetch_handle = Entrez.efetch(
                db="pubmed",
                rettype="medline",
                retmode="text",
                retstart=0,
                retmax=1000,
                webenv=search_results['WebEnv'],
                query_key=search_results['QueryKey']
            )
            
            # Process records
            papers = []
            for record in Entrez.parse(fetch_handle):
                paper_data = self._process_record(record)
                if paper_data:
                    papers.append(paper_data)
                    
            fetch_handle.close()
            
            self._logger.info(
                f"Found {len(papers)} papers with industry affiliations"
            )
            return papers
            
        except Exception as e:
            self._logger.error(f"Paper fetching error: {e}")
            raise

def setup_logging(debug: bool = False):
    """Configure logging with custom format and level"""
    logging.basicConfig(
        format='%(asctime)s [%(levelname)s] %(message)s',
        level=logging.DEBUG if debug else logging.INFO
    )

def main():
    """Command-line interface implementation"""
    parser = argparse.ArgumentParser(
        description='Analyze PubMed papers for industry affiliations'
    )
    parser.add_argument(
        'query',
        help='PubMed search query (supports full syntax)'
    )
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable detailed debug logging'
    )
    parser.add_argument(
        '-f', '--file',
        help='Output CSV filepath (prints to console if not specified)'
    )
    
    args = parser.parse_args()
    setup_logging(args.debug)
    logger = logging.getLogger(__name__)
    
    try:
        fetcher = PubMedResearchFetcher("your.email@example.com")
        papers = fetcher.fetch_papers(args.query)
        
        if not papers:
            logger.info("No relevant papers found")
            return
            
        # Convert to DataFrame
        results = pd.DataFrame([paper.to_dict() for paper in papers])
        
        # Handle output
        if args.file:
            output_path = Path(args.file)
            results.to_csv(output_path, index=False)
            logger.info(f"Results written to {output_path}")
        else:
            print(results.to_string())
            
    except Exception as e:
        logger.error(f"Execution failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
