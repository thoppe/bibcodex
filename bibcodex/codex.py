import pandas as pd
from typing import Dict
from .api import pubmed, semanticScholar, icite, doi2pmid


@pd.api.extensions.register_dataframe_accessor("bibcodex")
class Codex:

    # Reuse the components across all Codex instances
    pubmed = pubmed
    semanticScholar = semanticScholar
    icite = icite
    doi2pmid = doi2pmid

    def __init__(self, df):
        # Validation can happen here if needed
        self.df = df
        pass

    def clear(self) -> None:
        """
        Clears the cache for all databases.
        """
        API = [pubmed, semanticScholar, icite, doi2pmid]

        for api in API:
            api.clear()

    def validate(self) -> Dict:
        """
        Returns information about the number of valid PMIDs and DOIs
        """

        df = self.df
        n = len(df)

        if "pmid" in df:
            n_pmid = (~df["pmid"].isnull()).sum()
            n_pmid_missing = n - n_pmid
            n_pmid_invalid = n_pmid - self.valid_pmid_idx().sum()
        else:
            n_pmid = n_pmid_missing = n_pmid_invalid = 0

        if "doi" in df:
            n_doi = (~df["doi"].isnull()).sum()
            n_doi_missing = n - n_doi
            n_doi_invalid = n_doi - self.valid_doi_idx().sum()
        else:
            n_doi = n_doi_missing = n_doi_invalid = 0

        return {
            "n_rows": n,
            "n_pmid": n_pmid,
            "n_pmid_missing": n_pmid_missing,
            "n_pmid_invalid": n_pmid_invalid,
            "n_doi": n_doi,
            "n_doi_missing": n_doi_missing,
            "n_doi_invalid": n_doi_invalid,
        }

    def valid_doi_idx(self) -> pd.Series:
        """
        Returns a boolean series which marks invalid DOIs
        Missing values are considered invalid.
        """
        return self.df["doi"].apply(self.pubmed.check_doi)

    def valid_pmid_idx(self) -> pd.Series:
        """
        Returns a boolean series which marks invalid PMIDs
        Missing values are considered invalid.
        """
        return self.df["pmid"].apply(self.pubmed.check_pmid)

    def set_api_key(self, api: str, key: str) -> None:

        API = {
            "pubmed": pubmed,
            "semanticScholar": semanticScholar,
        }

        if api not in API:
            err = f"API key not implemented for {api}"
            raise NotImplementedError(err)

        API[api].api_key = key

    def download(self, api="pubmed", add_prefix=False, add_suffix=False):

        """
        Downloads (or pulls from the cache) data from an API (e.g. PubMed).
        Uses the current index (either a pmid or doi required).
        """

        method = self.df.index.name
        valid_methods = ["doi", "pmid"]
        df = self.df

        if df.index.name is None:
            msg = f"Index not set, should be one of {valid_methods}"
            raise TypeError(msg)

        if not pd.api.types.is_string_dtype(df.index):
            msg = f"Index of type is {method} is not a string (and should be!)"
            raise TypeError(msg)

        if method not in valid_methods:
            msg = f"Index must be one of {valid_methods} not {method}"
            raise TypeError(msg)

        atypes = ["pubmed", "icite", "semanticScholar", "doi2pmid"]
        if api not in atypes:
            raise NotImplementedError(f"API must be one of {atypes}")

        # Only collect data for the unique records that are not empty
        records = df.index.dropna().unique()

        if method == "doi":
            API = {
                "semanticScholar": self.semanticScholar,
                "doi2pmid": self.doi2pmid,
            }

        elif method == "pmid":
            records = [x for x in records if x != "nan"]

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

        # Drop fully missing values
        data = [row for row in data if row]

        # Cast the index to string
        data = pd.DataFrame(data)
        data[method] = data[method].astype(str)
        data = data.set_index(method)

        if add_prefix:
            data = data.add_prefix(f"{api}_")

        if add_suffix:
            data = data.add_suffix(f"_{api}")

        return data
