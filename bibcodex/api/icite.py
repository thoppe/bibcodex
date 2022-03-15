from .helpers import CachedDownloader, RemoteItemNotFound


class iCite_downloader(CachedDownloader):
    """
    Downloads and caches the information from iCite.
    Requires that the input be a non-negative integer.
    """

    name = "iCite"
    datatype = dict

    # has_multi_cache = True

    @CachedDownloader.cached
    def __call__(self, pmid: int) -> dict:

        self.validate_pmid(pmid)

        # if is_multi:
        #    return self.multi_download(pmid)

        url = "https://icite.od.nih.gov/api/pubs"
        r = self.download(url, {"pmids": pmid})

        assert r.ok

        js = r.json()

        assert "data" in js

        if len(js["data"]) < 1:
            raise RemoteItemNotFound(pmid, self.name)

        js = js["data"][0]
        return js

    """
    def multi_download(self, pmids):
        
        If a list of PMIDs are requested, return a dictionary for each PMID.
        iCite does not return a value if key is not found, so add None for missing values.
        

        # Require that all PMIDs are integers
        for k in pmids:
            if not isinstance(k, int):
                raise TypeError(f"{k} is expected to be an integer for {self.name}")

        str_pmids = ",".join(map(str, pmids))
        url = f"https://icite.od.nih.gov/api/pubs?pmids={str_pmids}"

        r = self.download(url)

        assert r.ok

        js = r.json()

        assert "data" in js

        data = {}
        for row in js["data"]:
            data[row["pmid"]] = row

        for k in pmids:
            if k not in data:
                data[k] = None

        return data
    """


downloader = iCite_downloader()
