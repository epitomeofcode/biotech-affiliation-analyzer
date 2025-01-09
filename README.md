# Research Papers Fetcher

A command-line tool to fetch research papers from PubMed with authors affiliated with pharmaceutical or biotech companies.

## Code Organization

The project is organized as follows:

```
get-papers-list/
├── src/
│   └── get_papers_list/
│       ├── __init__.py
│       └── main.py
├── tests/
│   └── test_main.py
├── pyproject.toml
└── README.md
```

- `src/get_papers_list/main.py`: Contains the main program logic
- `pyproject.toml`: Poetry configuration and dependency management
- `README.md`: Project documentation

## Installation

1. Make sure you have Poetry installed:
```bash
curl -sSL https://install.python-poetry.org | python3 -
```

2. Clone the repository:
```bash
git clone https://github.com/yourusername/get-papers-list.git
cd get-papers-list
```

3. Install dependencies:
```bash
poetry install
```

## Usage

The program can be run using the following command:

```bash
poetry run get-papers-list "your pubmed query"
```

### Command-line Options

- `-h, --help`: Display usage instructions
- `-d, --debug`: Enable debug logging
- `-f, --file FILE`: Save results to a CSV file instead of printing to console

### Example

```bash
poetry run get-papers-list "cancer therapy" -f results.csv
```

## Output Format

The program generates a CSV file with the following columns:

- PubmedID: Unique identifier for the paper
- Title: Title of the paper
- Publication Date: Date the paper was published
- Non-academic Author(s): Names of authors affiliated with companies
- Company Affiliation(s): Names of pharmaceutical/biotech companies
- Corresponding Author Email: Email address of the corresponding author

## Tools and Libraries Used

This project was built using the following tools and libraries:

- [Poetry](https://python-poetry.org/): Dependency management and packaging
- [Biopython](https://biopython.org/): Interface to PubMed/Entrez
- [Pandas](https://pandas.pydata.org/): Data manipulation and CSV handling
- [GitHub](https://github.com): Version control and hosting

The code structure and documentation were developed with assistance from GPT-4.

## License

MIT License
