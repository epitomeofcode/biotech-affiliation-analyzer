# PubMed Industry Collaboration Analyzer

An advanced tool for analyzing research papers from PubMed to identify industry-academic collaborations and pharmaceutical/biotech company affiliations.

## Features

- Advanced affiliation analysis using multi-factor classification
- Robust date parsing for various PubMed formats
- Structured data handling using dataclasses
- Comprehensive logging system
- Type-annotated codebase
- Configurable output formats

## Project Structure

```
pubmed-industry-analyzer/
├── src/
│   └── pubmed_industry_analyzer/
│       ├── __init__.py
│       ├── main.py
│       ├── affiliation.py
│       └── parser.py
├── tests/
│   ├── __init__.py
│   ├── test_affiliation.py
│   └── test_parser.py
├── pyproject.toml
└── README.md
```

## Installation

1. Ensure Poetry is installed on your system:
```bash
pipx install poetry
```

2. Clone and set up the project:
```bash
git clone https://github.com/yourusername/pubmed-industry-analyzer.git
cd pubmed-industry-analyzer
poetry install
```

## Usage

Basic usage:
```bash
poetry run analyze-pubmed "your search query"
```

Advanced options:
```bash
poetry run analyze-pubmed "cancer immunotherapy" -f results.csv -d
```

### Command Arguments

- Required: PubMed search query (supports full PubMed syntax)
- `-d/--debug`: Enable detailed debug logging
- `-f/--file PATH`: Save results to specified CSV file

## Output Format

The tool generates a CSV with the following columns:

- `PubmedID`: Unique identifier
- `Title`: Paper title
- `Publication_Date`: Standardized publication date
- `Non_academic_Authors`: Industry-affiliated authors
- `Company_Affiliations`: Identified company names
- `Corresponding_Author_Email`: Contact information

## Development Notes

This project implements several unique features:

1. Multi-factor affiliation analysis
   - Company name pattern matching
   - Industry keyword scoring
   - Academic institution filtering

2. Robust data structures
   - Type-safe dataclasses
   - Structured error handling
   - Clear separation of concerns

3. Quality assurance
   - Type checking with mypy
   - Style enforcement with black
   - Linting with pylint

## Dependencies

Core libraries:
- Biopython: PubMed API integration
- Pandas: Data manipulation
- Python-dateutil: Date parsing

Development tools:
- Poetry: Dependency management
- PyTest: Testing framework
- MyPy: Static type checking

## Contributing

1. Fork the repository
2. Create a feature branch
3. Implement your changes
4. Add tests for new functionality
5. Submit a pull request

## License

MIT License - See LICENSE file for details
