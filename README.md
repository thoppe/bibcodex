# bibcodex
_🚧 WORK IN PROGRESS 🚧_

Library to access, analyze, and display bibliographic information. 

## Examples

Import the `pandas` and `bibcodex` together and load a dataframe:
```python
import bibcodex
import pandas as pd

# You should always cast your search variables (pmid, doi) to str.
df = pd.read_csv("data/sample_data.csv", dtype={'pmid':str})
```

Valid download methods are: `icite`, `doi2pmid`, `semanticScholar`, or `pubmed`:

```python
# Set the index to search query
df = df.set_index("doi")

# Download the information, and combine it with the original dataframe:
info = df.codex.download('semanticScholar')
print(df.combine_first(info[["title"]]))

"""
doi                      title                                                           
10.1001/jama.2017.18444  Progressive Massive Fibrosis in Coal Miners Fr...
10.1001/jama.2018.0126   Birth Defects Potentially Related to Zika Viru...
10.1001/jama.2018.0708   Association Between Estimated Cumulative Vacci...
10.1001/jama.2018.10488  Electronic Cigarette Sales in the United State...
"""
```

All search queries are cached locally in `./cache`. To clear the cache use:

```python
df.codex.clear()
```


| API  | Returned fields |
| ------------- | ------------- |
| `pubmed`           | title issue pages abstract journal authors pubdate mesh_terms publication_types chemical_list keywords doi references delete languages vernacular_title affiliations pmc other_id medline_ta nlm_unique_id issn_linking country  |
| `semanticScholar`  | abstract arxivId authors citationVelocity citations corpusId fieldsOfStudy influentialCitationCount isOpenAccess isPublisherLicensed is_open_access is_publisher_licensed numCitedBy numCiting paperId references s2FieldsOfStudy title topics url venue year  |
| `icite`            | year title authors journal is_research_article relative_citation_ratio nih_percentile human animal molecular_cellular apt is_clinical citation_count citations_per_year expected_citations_per_year field_citation_rate provisional x_coord y_coord cited_by_clin cited_by references doi  |
| `doi2pmid`         | live status errmsg pmcid pmid versions  |




## Roadmap

- [x] API access: Pubmed (Parsed MEDLINE data)
- [x] API access: Semantic Scholar (PMID)
- [x] API access: iCite
- [x] API access: Semantic Scholar (DOI)
- [x] API access: DOI to PMID NLM www.ncbi.nlm.nih.gov/pmc/tools/idconv/
- [ ] API access: Pubmed (XML)
- [ ] API access: arXiv
- [ ] API access: CoLIL
- [x] API access, validation of input
- [x] API access, multi item requests
- [x] API access, chunking
- [ ] API access, include status_code in download results 
- [ ] API access, better error handling
- [ ] API caching, clearing
- [ ] API caching, in-class domain (what did I mean here??)
- [x] Codex, validate PMID
- [x] Codex, validate DOI
- [x] Codex, build dataframe from items
- [ ] Testing harness
- [x] CI linting
- [ ] pypi library (poetry?)
- [x] README with examples
- [ ] Embedding functions (SPECTER)
- [ ] Clustering
- [ ] Visualization (streamlit)

