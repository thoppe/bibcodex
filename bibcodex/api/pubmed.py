from Bio import Entrez
from .helpers import CachedDownloader, cached, RemoteItemNotFound


class PubMed_downloader(CachedDownloader):

    name = "pubmed"
    datatype = str
    api_key = ""

    def validate(self, pmid):
        if not isinstance(pmid, int):
            raise TypeError(f"{self.name} expected an integer, called with {pmid}")

    def download(self, pmid):
        # Need custom download, since this uses Bio.Entrez
        self.validate(pmid)

        # Fill this in for PubMed team
        Entrez.email = "https://github.com/thoppe/bibcodex"

        # If this is filled, we can go faster
        Entrez.api_key = PubMed_downloader.api_key

        res = Entrez.efetch(db="pubmed", id=[str(pmid)], retmode="xml")
        return res.read()

    @cached
    def __call__(self, pmid):
        r = self.download(pmid)
        xml = r.decode()

        # Entrez returns an empty-ish XML that looks like this
        if "<PubmedArticleSet></PubmedArticleSet>" in xml:
            raise RemoteItemNotFound(pmid, self.name)

        return xml


downloader = PubMed_downloader()
