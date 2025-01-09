import argparse
import sys
import logging
from typing import List, Dict, Optional
from Bio import Entrez
import pandas as pd
import csv
from pathlib import Path

# Configure logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger = logging.getLogger(__name__)

class PapersFetcher:
    def __init__(self, email: str):
        """Initialize the PapersFetcher with user email for PubMed access."""
        Entrez.email = email
        self.company_keywords = [
            'pharma', 'biotech', 'therapeutics', 'pharmaceuticals', 
            'biosciences', 'laboratories', 'medicine'
        ]
        self.academic_keywords = [
            'university', 'college', 'institute', 'school',
            'hospital', 'medical center', 'clinic'
        ]

    def is_company_affiliation(self, affiliation: str) -> bool:
        """Determine if an affiliation is from a pharmaceutical/biotech company."""
        if not affiliation:
            return False
            
        affiliation_lower = affiliation.lower()
        
        # Check for company indicators
        has_company_keyword = any(keyword in affiliation_lower 
                                for keyword in self.company_keywords)
                                
        # Check for academic indicators
        is_academic = any(keyword in affiliation_lower 
                         for keyword in self.academic_keywords)
                         
        return has_company_keyword and not is_academic

    def fetch_papers(self, query: str) -> List[Dict]:
        """Fetch papers from PubMed based on the query."""
        logger.info(f"Searching PubMed with query: {query}")
        
        try:
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
            record = Entrez.read(handle)
            handle.close()

            if not record['IdList']:
                logger.warning("No papers found matching the query")
                return []

            # Fetch paper details
            handle = Entrez.efetch(db="pubmed", id=record['IdList'], 
                                 rettype="medline", retmode="text")
            records = Entrez.parse(handle)
            
            papers = []
            for record in records:
                paper = self._process_paper(record)
                if paper:
                    papers.append(paper)
                    
            logger.info(f"Found {len(papers)} papers with company affiliations")
            return papers
            
        except Exception as e:
            logger.error(f"Error fetching papers: {e}")
            raise

    def _process_paper(self, record: Dict) -> Optional[Dict]:
        """Process a single paper record and extract relevant information."""
        try:
            # Extract authors and affiliations
            authors = record.get('AU', [])
            affiliations = record.get('AD', [])
            
            # Process affiliations
            company_authors = []
            company_affiliations = set()
            
            for author, affiliation in zip(authors, affiliations):
                if self.is_company_affiliation(affiliation):
                    company_authors.append(author)
                    company_affiliations.add(affiliation)
            
            if not company_authors:
                return None
                
            # Extract corresponding author email
            email = None
            for field in record.get('FAU', []):
                if 'corresponding author' in field.lower():
                    email_match = re.search(r'[\w\.-]+@[\w\.-]+', field)
                    if email_match:
                        email = email_match.group(0)
                        break
            
            return {
                'PubmedID': record.get('PMID', ''),
                'Title': record.get('TI', ''),
                'Publication_Date': record.get('DP', ''),
                'Non_academic_Authors': '; '.join(company_authors),
                'Company_Affiliations': '; '.join(company_affiliations),
                'Corresponding_Author_Email': email or ''
            }
            
        except Exception as e:
            logger.error(f"Error processing paper {record.get('PMID', '')}: {e}")
            return None

def main():
    """Main entry point for the command-line program."""
    parser = argparse.ArgumentParser(
        description='Fetch research papers from PubMed with company affiliations'
    )
    parser.add_argument('query', help='PubMed search query')
    parser.add_argument('-d', '--debug', action='store_true', 
                       help='Enable debug logging')
    parser.add_argument('-f', '--file', help='Output file path (CSV)')
    args = parser.parse_args()

    # Configure debug logging if requested
    if args.debug:
        logger.setLevel(logging.DEBUG)

    try:
        # Initialize fetcher with your email
        fetcher = PapersFetcher("your.email@example.com")
        papers = fetcher.fetch_papers(args.query)

        if not papers:
            logger.info("No papers found with company affiliations")
            return

        # Create DataFrame
        df = pd.DataFrame(papers)

        # Output handling
        if args.file:
            output_path = Path(args.file)
            df.to_csv(output_path, index=False)
            logger.info(f"Results saved to {output_path}")
        else:
            print(df.to_string())

    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
