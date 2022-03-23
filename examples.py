import pandas as pd
from bibcodex.api import pubmed, semanticScholar, icite


p1 = 27034526
p0 = 31633016

# iCite returns a dictionary of information
print(icite(p0))

# {'pmid': 31633016, 'year': 2019, 'title': 'Topic choice contributes to the lower rate of NIH awards to African-American/black scientists.', 'authors': 'Travis A Hoppe, Aviva Litovitz, Kristine A Willis, Rebecca A Meseroll, Matthew J Perkins, B Ian Hutchins, Alison F Davis, Michael S Lauer, Hannah A Valantine, James M Anderson, George M Santangelo', 'journal': 'Sci Adv', 'is_research_article': 'Yes', 'relative_citation_ratio': 19.41, 'nih_percentile': 99.5, 'human': 1.0, 'animal': 0.0, 'molecular_cellular': 0.0, 'apt': 0.95, 'is_clinical': 'No', 'citation_count': 172, 'citations_per_year': 57.333333333333336,

# PubMed returns the raw XML. Multiple PMIDs can be passed at once
print(pubmed([p0, p1]))

# ['<PubmedArticle><MedlineCitation Owner="NLM" Status="MEDLINE"><PMID Version="1">31633016</PMID><DateCompleted><Year>2020</Year><Month>05</Month><Day>19</Day></DateCompleted><DateRevised><Year>2021</Year><Month>01</Month><Day>10</Day></DateRevised><Article PubModel="Electronic-eCollection"><Journal><ISSN IssnType="Electronic">2375-2548</ISSN><JournalIssue CitedMediu

# Semantic Scholar returns *very* detialed information, but API has limited requests
# bibcodex caches previous results so only one API call is needed
print(semanticScholar(p0))

# [{'abstract': 'Topic choice is a previously unappreciated contributor to the lower rate of NIH awards to AA/B scientists. Despite efforts to promote diversity in the biomedical workforce, there remains a lower rate of funding of National Institutes of Health R01 applications submitted by African-American/black (AA/B) scientists relative to white scientists. To identify underlying causes of this funding gap, we analyzed six stages of the application process from 2011 to 2015 and found that disparate outcomes arise at three of the six: decision to discuss, impact score assignment, and a previously unstudied stage, topic choice. Notably, AA/B applicants tend to propose research on topics with lower award rates. These topics include research at the community and population level, as opposed to more fundamental and mechanistic investigations; the latter tend to have higher award rates. Topic choice alone accounts for over 20% of the funding gap after controlling for multiple variables, including the applicant’s prior achievements. Our findings can be used to inform interventions designed to close the funding gap.', 'arxivId': None, 'authors': [{'authorId': '47000911', 'name': 'Travis Hoppe', 'url': 'https://www.semanticscholar.org/author/47000911'}, {'authorId': '5518590',
