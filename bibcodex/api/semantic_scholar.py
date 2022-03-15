from .helpers import CachedDownloader, RemoteItemNotFound


class SemanticScholar_downloader(CachedDownloader):
    """
    Downloads and caches the output from api.semanticscholar.org
    Requires that the input be a non-negative integer.
    """

    name = "SemanticScholar"
    datatype = dict

    @CachedDownloader.cached
    def __call__(self, pmid: int) -> dict:

        self.validate_pmid(pmid)

        # Custom URL
        url = (
            f"http://s2-public-api-prod.us-west-2.elasticbeanstalk.com/"
            f"v1/paper/"
            f"PMID:{pmid}?include_unknown_references=true"
        )

        # Public facing API
        url = f"https://api.semanticscholar.org/v1/paper/PMID:{pmid}"

        r = self.download(url)

        if r.status_code in [
            404,
        ]:
            raise RemoteItemNotFound(pmid, self.name)

        js = r.json()

        return js


downloader = SemanticScholar_downloader()
