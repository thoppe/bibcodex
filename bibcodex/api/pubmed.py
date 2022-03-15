from Bio import Entrez
from typing import List, Union, Dict

from .helpers import CachedDownloader, RemoteItemNotFound

import pandas as pd
from bs4 import BeautifulSoup


class PubMed_downloader(CachedDownloader):
    """
    Downloads and caches the raw XML from PubMed.
    Requires that the input be a non-negative integer.
    """

    name = "pubmed"
    datatype = str

    api_key = ""

    def download(self, pmids: int) -> str:
        # Need custom download, since this uses Bio.Entrez

        # Fill this in for PubMed team
        Entrez.email = "https://github.com/thoppe/bibcodex"

        # If this is filled, we can go faster
        Entrez.api_key = PubMed_downloader.api_key

        # Call the Entrez downloader
        res = Entrez.efetch(db="pubmed", id=pmids, retmode="xml")

        return res.read()

    @CachedDownloader.cached
    def __call__(self, pmids: Union[int, List]) -> Dict[str, datatype]:

        # Validate the input datatypes
        [self.validate_pmid(p) for p in pmids]

        # Download the raw data
        xml = self.download(pmids)

        # Parse the returned XML and extract each result
        soup = BeautifulSoup(xml, "xml")

        data = {}
        for block in soup.find_all("PubmedArticle"):
            info = block.PubmedData.ArticleIdList
            pmid = info.find(attrs={"IdType": "pubmed"})
            pmid = pmid.get_text().strip()
            data[pmid] = str(block)

        return data


downloader = PubMed_downloader()
