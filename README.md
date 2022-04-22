# bibcodex

Library to access, analyze, and display bibliographic information.
__ 🚧 WORK IN PROGRESS 🚧 __

## Examples

Import the `pandas` and `bibcodex` together:
```python
import bibcodex
import pandas as pd

df = pd.read_csv("data/sample_data.csv", dtype={'pmid':str})
```

You should always cast your search variables (pmid, doi) to str.

```python
df = df.set_index("doi")

info = df.codex.download("doi2pmid")
print(df.combine_first(info[["pmid"]]))
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
- [ ] README with examples
- [ ] Embedding functions (SPECTER)
- [ ] Clustering
- [ ] Visualization (streamlit)

