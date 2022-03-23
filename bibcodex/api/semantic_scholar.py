from typing import List, Dict
from .helpers import CachedDownloader


class SemanticScholar_downloader(CachedDownloader):
    """
    Downloads and caches the output from api.semanticscholar.org
    Requires that the input be a non-negative integer.

    API reference:
    https://www.semanticscholar.org/product/api#Fetch-Paper
    """

    name = "semanticScholar"
    datatype = dict
    chunksize = 1

    def access_API(self, records: List[int], url: str) -> Dict[str, datatype]:

        if len(records) > 1:
            raise ValueError("Semantic Scholar can not do multi downloads")

        record = str(records[0])
        headers = {}

        # Apply the API if it exists
        if self.api_key:
            headers["x-api-key"] = self.api_key

        r = self.download(url + record, headers=headers)

        if r.status_code in [404]:
            return {}

        return r.json(), record

    @CachedDownloader.cached
    def get_from_PMIDs(self, pmids: List[int]) -> Dict[str, datatype]:

        url = "https://api.semanticscholar.org/v1/paper/PMID:"
        results, pmid = self.access_API(pmids, url)

        # Add in the pmid as Semantic Scholar doesn't return it
        results["pmid"] = pmid

        return {str(pmid): results}

    @CachedDownloader.cached
    def get_from_DOIs(self, dois: List[int]) -> Dict[str, datatype]:

        url = "https://api.semanticscholar.org/v1/paper/"
        results, doi = self.access_API(dois, url)

        return {str(doi): results}


downloader = SemanticScholar_downloader()
