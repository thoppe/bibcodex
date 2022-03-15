from typing import List, Union, Dict
from .helpers import CachedDownloader, RemoteItemNotFound


class SemanticScholar_downloader(CachedDownloader):
    """
    Downloads and caches the output from api.semanticscholar.org
    Requires that the input be a non-negative integer.
    """

    name = "SemanticScholar"
    datatype = dict
    chunksize = 1

    @CachedDownloader.cached
    def __call__(self, pmids: Union[int, List]) -> Dict[str, datatype]:

        # Validate the input datatypes
        [self.validate_pmid(p) for p in pmids]

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
