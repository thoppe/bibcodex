from typing import List, Dict
from .helpers import CachedDownloader


# Example API
# https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my_email@example.com&ids=PMC3531190,PMC3245039


class PMC_IdConverter_downloader(CachedDownloader):
    """
    Downloads and caches output from ID Converter converting DOIs to PMIDs.
    https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

    """

    name = "idConverter"
    datatype = dict
    chunksize = 100

    @CachedDownloader.cached
    def get_from_DOIs(self, dois: List[int]) -> Dict[str, datatype]:

        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        params = {
            "email": "travis.hoppe+bibcodex@gmail.com",
            "tool": "https://github.com/thoppe/bibcodex",
            "idtype": "doi",
            "format": "json",
            "ids": ",".join(dois),
        }

        result = self.download(url, params=params)

        data = {}
        for row in result.json()["records"]:
            data[row["doi"]] = row

        return data


downloader = PMC_IdConverter_downloader()
