from Bio import Entrez
from .helpers import CachedDownloader, cached


class PubMed_downloader(CachedDownloader):

    name = "pubmed"
    datatype = str
    api_key = ""

    def download(self, pmid):
        # Need custom download, since this uses Bio.Entrez

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
        return xml


downloader = PubMed_downloader()
