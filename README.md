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

Valid download methods are: `icite`, `doi2pmid`, 'semanticScholar`, or `pubmed`:

```python
# Set the index to search query
df = df.set_index("doi")

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

