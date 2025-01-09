"""
PubMed Industry Affiliation Analyzer
A tested implementation for analyzing research papers with industry affiliations.
"""
import argparse
import sys
import logging
import time
from typing import List, Dict, Optional, Set
from dataclasses import dataclass
from datetime import datetime
from Bio import Entrez, Medline
import pandas as pd
import re
from pathlib import Path

@dataclass
class PaperAuthor:
    """Author information container"""
    name: str
    affiliation: str
    email: Optional[str] = None
    is_corresponding: bool = False

@dataclass
class Paper:
    """Paper information container"""
    pubmed_id: str
    title: str
    pub_date: str
    authors: List[PaperAuthor]

    def to_dict(self) -> Dict:
        """Convert to dictionary format for CSV export"""
        industry_authors = [a for a in self.authors 
                          if a.affiliation]  # Only include authors with affiliations
        return {
            'PubmedID': self.pubmed_id,
            'Title': self.title,
            'Publication_Date': self.pub_date,
            'Non_academic_Authors': '; '.join(a.name for a in industry_authors),
            'Company_Affiliations': '; '.join(set(a.affiliation for a in industry_authors)),
            'Corresponding_Author_Email': next((a.email for a in self.authors if a.is_corresponding), '')
        }

class AffiliationChecker:
    """Analyzes author affiliations"""
    
    def __init__(self):
        self.company_indicators: Set[str] = {'inc', 'corp', 'ltd', 'llc', 'gmbh', 'co', 'pharmaceuticals', 'biotech', 'therapeutics'}
        self.academic_indicators: Set[str] = {'university', 'college', 'institute', 'school','hospital', 'clinic', 'medical center' }

    def is_company(self, affiliation: str) -> Optional[str]:
        """
        Check if affiliation is a company
        Returns cleaned company name if true, None otherwise
        """
        if not affiliation:
            return None

        affiliation_lower = affiliation.lower()
        
        # Quick reject for academic institutions
        if any(ind in affiliation_lower for ind in self.academic_indicators):
            return None

        # Check for company indicators
        if any(ind in affiliation_lower for ind in self.company_indicators):
            # Extract company name - take the first part before common separators
            parts = re.split(r'[,;]', affiliation)
            return parts[0].strip()

        return None

class PubMedFetcher:
    """Handles PubMed API interactions"""
    
    def __init__(self, email: str):
        if not email or '@' not in email:
            raise ValueError("Valid email required for PubMed API access")
        self.email = email
        Entrez.email = email
        self.affiliation_checker = AffiliationChecker()
        self.logger = logging.getLogger(__name__)

    def _extract_email(self, text: str) -> Optional[str]:
        """Extract email from text"""
        if not text:
            return None
        matches = re.findall(r'[\w\.-]+@[\w\.-]+\.\w+', text)
        return matches[0] if matches else None

    def _parse_date(self, medline_date: str) -> str:
        """Parse PubMed date format"""
        try:
            # Handle various date formats
            if not medline_date:
                return ""
                
            # Try to parse YYYY Mon DD format
            match = re.match(r'(\d{4}) (\w+)( \d+)?', medline_date)
            if match:
                year = match.group(1)
                month = match.group(2)
                day = match.group(3).strip() if match.group(3) else "1"
                try:
                    date_obj = datetime.strptime(f"{year} {month} {day}", "%Y %b %d")
                    return date_obj.strftime("%Y-%m-%d")
                except ValueError:
                    return year
            return medline_date
        except Exception as e:
            self.logger.warning(f"Date parsing error: {e}")
            return medline_date

    def _process_record(self, record: Medline.Record) -> Optional[Paper]:
        """Process a single PubMed record"""
        try:
            # Basic paper info
            pubmed_id = record.get('PMID', '')
            title = record.get('TI', '').replace('\n', ' ')
            
            # Handle authors and affiliations
            authors = []
            affiliations = record.get('AD', '').split(';')
            author_names = record.get('FAU', [])
            
            # Match authors with affiliations
            for idx, name in enumerate(author_names):
                affiliation = affiliations[idx] if idx < len(affiliations) else ''
                company_name = self.affiliation_checker.is_company(affiliation)
                
                if company_name:
                    is_corresponding = 'correspond' in affiliation.lower()
                    email = self._extract_email(affiliation) if is_corresponding else None
                    
                    authors.append(PaperAuthor(
                        name=name,
                        affiliation=company_name,
                        email=email,
                        is_corresponding=is_corresponding
                    ))
            
            if not authors:  # Skip if no industry authors found
                return None
                
            return Paper(
                pubmed_id=pubmed_id,
                title=title,
                pub_date=self._parse_date(record.get('DP', '')),
                authors=authors
            )
            
        except Exception as e:
            self.logger.error(f"Error processing record {record.get('PMID', '')}: {e}")
            return None

    def fetch_papers(self, query: str) -> List[Paper]:
        """Fetch papers matching query"""
        self.logger.info(f"Searching PubMed: {query}")
        papers = []
        
        try:
            # Initial search
            handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
            results = Entrez.read(handle)
            handle.close()
            
            if not results['IdList']:
                self.logger.info("No results found")
                return papers
                
            # Fetch details with rate limiting
            for paper_id in results['IdList']:
                time.sleep(0.34)  # Rate limit: max 3 requests per second
                try:
                    fetch_handle = Entrez.efetch(
                        db="pubmed",
                        id=paper_id,
                        rettype="medline",
                        retmode="text"
                    )
                    record = Medline.read(fetch_handle)
                    fetch_handle.close()
                    
                    paper = self._process_record(record)
                    if paper:
                        papers.append(paper)
                        
                except Exception as e:
                    self.logger.error(f"Error fetching paper {paper_id}: {e}")
                    continue
                    
            self.logger.info(f"Found {len(papers)} papers with industry affiliations")
            return papers
            
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            raise

def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze PubMed papers for industry affiliations'
    )
    parser.add_argument('query', help='PubMed search query')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('-f', '--file', help='Output CSV file path')
    parser.add_argument('-e', '--email', required=True, help='Email for PubMed API access')
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    logger = logging.getLogger(__name__)
    
    try:
        fetcher = PubMedFetcher(args.email)
        papers = fetcher.fetch_papers(args.query)
        
        if not papers:
            logger.info("No matching papers found")
            return
            
        # Convert to DataFrame
        df = pd.DataFrame([paper.to_dict() for paper in papers])
        
        # Handle output
        if args.file:
            output_path = Path(args.file)
            df.to_csv(output_path, index=False)
            logger.info(f"Results saved to {output_path}")
        else:
            print(df.to_string())
            
    except Exception as e:
        logger.error(f"Program failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
