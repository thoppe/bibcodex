from typing import List, Dict
from .helpers import CachedDownloader


class iCite_downloader(CachedDownloader):
    """
    Downloads and caches the information from iCite.
    Requires that the input be a non-negative integer.
    """

    name = "iCite"
    datatype = dict
    chunksize = 100

    @CachedDownloader.cached
    def get_from_PMIDs(self, pmids: List[int]) -> Dict[str, datatype]:
        str_pmids = ",".join(map(str, pmids))
        url = f"https://icite.od.nih.gov/api/pubs?pmids={str_pmids}"

        r = self.download(url)
        assert r.ok

        js = r.json()

        data = {}
        for item in js["data"]:
            pmid = str(item["pmid"])
            data[pmid] = item

        return data


downloader = iCite_downloader()
