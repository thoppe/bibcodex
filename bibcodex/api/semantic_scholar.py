from typing import List, Dict
from .helpers import CachedDownloader


class SemanticScholar_downloader(CachedDownloader):
    """
    Downloads and caches the output from api.semanticscholar.org
    Requires that the input be a non-negative integer.
    """

    name = "SemanticScholar"
    datatype = dict
    chunksize = 1

    @CachedDownloader.cached
    def get_from_PMIDs(self, pmids: List[int]) -> Dict[str, datatype]:

        if len(pmids) > 1:
            raise ValueError("Semantic Scholar can not do multi downloads")

        pmid = pmids[0]

        # Public facing API
        url = f"https://api.semanticscholar.org/v1/paper/PMID:{pmid}"

        r = self.download(url)

        if r.status_code in [404]:
            return {}

        return {str(pmid): r.json()}


downloader = SemanticScholar_downloader()
