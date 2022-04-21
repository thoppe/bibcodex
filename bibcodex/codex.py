import pandas as pd
from typing import Dict
from .api import pubmed, semanticScholar, icite, idConverter


class Codex(pd.DataFrame):

    # Reuse the components across all Codex instances
    pubmed = pubmed
    semanticScholar = semanticScholar
    icite = icite
    idConverter = idConverter

    def validate(self) -> Dict:
        """
        Returns information about the number of valid PMIDs and DOIs
        """
        if "pmid" in self:
            n_pmid = sum(self["pmid"].notnull())
            n_pmid_missing = len(self) - n_pmid
            n_pmid_invalid = sum(self.invalid_pmid)

        else:
            n_pmid = n_pmid_invalid = n_pmid_invalid = 0

        if "doi" in self:
            n_doi = sum(self["doi"].notnull())
            n_doi_missing = len(self) - n_doi
            n_doi_invalid = sum(self.invalid_doi)
        else:
            n_doi = n_doi_missing = n_doi_invalid = 0

        return {
            "n_rows": len(self),
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

    def clear_cache(self) -> None:
        """
        Clears the cache for all databases.
        """
        API = [self.pubmed, self.semanticScholar, self.icite, self.idConverter]
        for api in API:
            api.clear()

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
