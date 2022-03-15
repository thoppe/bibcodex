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

    def download(self, pmid: int) -> str:
        # Need custom download, since this uses Bio.Entrez

        # Fill this in for PubMed team
        Entrez.email = "https://github.com/thoppe/bibcodex"

        # If this is filled, we can go faster
        Entrez.api_key = PubMed_downloader.api_key

        res = Entrez.efetch(db="pubmed", id=[str(pmid)], retmode="xml")

        return res.read()

    @CachedDownloader.cached
    def __call__(self, pmid) -> datatype:

        self.validate_pmid(pmid)

        xml = self.get_from_pmid(pmid)

        # Entrez returns an empty-ish XML that looks like this
        if "<PubmedArticleSet></PubmedArticleSet>" in xml:
            raise RemoteItemNotFound(pmid, self.name)

        return xml

    def get_from_pmid(self, pmid: int) -> datatype:
        r = self.download(pmid)
        xml = r.decode()
        return xml


downloader = PubMed_downloader()
