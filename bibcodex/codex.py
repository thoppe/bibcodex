import pandas as pd
from typing import Dict
from .api import pubmed, semanticScholar, icite, idConverter


@pd.api.extensions.register_dataframe_accessor("codex")
class Codex:

    # Reuse the components across all Codex instances
    pubmed = pubmed
    semanticScholar = semanticScholar
    icite = icite
    idConverter = idConverter

    def __init__(self, df):
        # Validation can happen here if needed
        self.df = df
        pass

    def clear(self) -> None:
        """
        Clears the cache for all databases.
        """
        API = [self.pubmed, self.semanticScholar, self.icite, self.idConverter]

        for api in API:
            api.clear()

    def validate(self) -> Dict:
        """
        Returns information about the number of valid PMIDs and DOIs
        """
        n = len(self.df)

        if "pmid" in self.df:
            invalid_pmid = self.df["pmid"].notnull()
            n_pmid = sum(~invalid_pmid)
            n_pmid_missing = n - n_pmid
            n_pmid_invalid = sum(invalid_pmid)

        else:
            n_pmid = n_pmid_missing = n_pmid_invalid = 0

        if "doi" in self.df:
            invalid_doi = self.df["doi"].notnull()
            n_doi = sum(~invalid_doi)
            n_doi_missing = n - n_doi
            n_doi_invalid = sum(invalid_doi)
        else:
            n_doi = n_doi_missing = n_doi_invalid = 0

        return {
            "n_rows": n,
            "n_pmid": n_pmid,
            "n_pmid_missing": n_pmid_missing,
            "n_pmid_invalid": n_pmid_invalid - n_pmid_missing,
            "n_doi": n_doi,
            "n_doi_missing": n_doi_missing,
            "n_doi_invalid": n_doi_invalid - n_doi_missing,
        }

    @property
    def invalid_doi(self) -> pd.Series:
        """
        Returns a boolean series which marks invalid DOIs
        Missing values are considered invalid.
        """
        return ~self["doi"].apply(self.pubmed.check_doi)

    @property
    def invalid_pmid(self) -> pd.Series:
        """
        Returns a boolean series which marks invalid PMIDs
        Missing values are considered invalid.
        """
        return ~self["pmid"].apply(self.pubmed.check_pmid)

    def set_api_key(self, api: str, key: str) -> None:

        API = {
            "pubmed": self.pubmed,
            "semanticScholar": self.semanticScholar,
        }

        if api not in API:
            err = f"API key not implemented for {api}"
            raise NotImplementedError(err)

        API[api].api_key = key

    def enrich(
        self, method="pmid", api="pubmed", add_prefix=True, add_suffix=False
    ):

        """
        Downloads (or pulls from the cache) data from an API (e.g. PubMed)
        and returns a Codex dataframe.
        """

        etypes = ["pmid", "doi"]
        if method not in etypes:
            raise NotImplementedError(f"enrich method must be one of {etypes}")

        atypes = ["pubmed", "icite", "semanticScholar", "idConverter"]
        if api not in atypes:
            raise NotImplementedError(f"enrich API must be one of {atypes}")

        if method not in self:
            raise ValueError(f"Method column {method} not in codex")

        # Only collect data for the records that are not empty
        records = self[method].dropna()

        if method == "doi":
            API = {
                "semanticScholar": self.semanticScholar,
                "idConverter": self.idConverter,
            }

        elif method == "pmid":
            # Convert method column to string representation or None
            self[method] = [
                str(int(x)) if not pd.isnull(x) else None for x in self[method]
            ]

            # Only collect data for the records that are not empty
            records = self[method].dropna()

            API = {
                "pubmed": self.pubmed,
                "icite": self.icite,
                "semanticScholar": self.semanticScholar,
            }

        # Check if the API has the method implemented
        if api not in API:
            err = f"{method} not implemented for {api}"
            raise NotImplementedError(err)

        # Call the API
        data = API[api](records, method=method)

        # Cast the index to string
        data = pd.DataFrame(data)
        data[method] = data[method].astype(str)
        data = data.set_index(method)

        if add_prefix:
            data = data.add_prefix(f"{api}_")

        if add_suffix:
            data = data.add_suffix(f"_{api}")

        # Merge the results together
        df = self.set_index(method).combine_first(data).reset_index()

        # Make sure we return a codex object
        self = Codex(df)

        # Return self to help chaining
        return self


"""
def read_csv(*args, **kwargs) -> Codex:
    # Wrap read_csv so Codex can call it directly
    return Codex(pd.read_csv(*args, **kwargs))
"""
