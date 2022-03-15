from Bio import Entrez
from .helpers import CachedDownloader, RemoteItemNotFound


class PubMed_downloader(CachedDownloader):
    """
    Downloads and caches the raw XML from PubMed.
    Requires that the input be a non-negative integer.
    """

    name = "pubmed"
    datatype = str
    api_key = ""

    def validate(self, pmid: int):

        if not isinstance(pmid, int) or pmid <= 0:
            err = f"{self.name} expected an non-negative integer (PMID), called with {pmid}"
            raise TypeError(err)

    def download(self, pmid: int) -> str:
        # Need custom download, since this uses Bio.Entrez
        self.validate(pmid)

        # Fill this in for PubMed team
        Entrez.email = "https://github.com/thoppe/bibcodex"

        # If this is filled, we can go faster
        Entrez.api_key = PubMed_downloader.api_key

        res = Entrez.efetch(db="pubmed", id=[str(pmid)], retmode="xml")

        return res.read()

    @CachedDownloader.cached
    def __call__(self, pmid) -> str:
        r = self.download(pmid)
        xml = r.decode()

        # Entrez returns an empty-ish XML that looks like this
        if "<PubmedArticleSet></PubmedArticleSet>" in xml:
            raise RemoteItemNotFound(pmid, self.name)

        return xml


downloader = PubMed_downloader()
