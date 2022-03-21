import pandas as pd
from .api import pubmed, semantic_scholar, icite


class Codex(pd.DataFrame):
    def enrich(
        self, method="pmid", api="pubmed", add_prefix=True, add_suffix=False
    ):

        etypes = ["pmid", "doi"]
        if method not in etypes:
            raise NotImplementedError(f"enrich method must be one of {etypes}")

        atypes = ["pubmed", "icite", "semantic_scholar"]
        if api not in atypes:
            raise NotImplementedError(f"enrich API must be one of {atypes}")

        if method not in self:
            raise ValueError(f"Method column {method} not in codex")

        # Only collect data for the records that are not empty
        records = self[method].dropna()

        if method == "doi":
            methods = {
                "semantic_scholar": semantic_scholar,
            }

        elif method == "pmid":
            # Convert method column to string representation or None
            self[method] = [
                str(int(x)) if not pd.isnull(x) else None for x in self[method]
            ]

            # Only collect data for the records that are not empty
            records = self[method].dropna()

            methods = {
                "pubmed": pubmed,
                "icite": icite,
                "semantic_scholar": semantic_scholar,
            }

        # Check if the API has the method implemented
        if api not in methods:
            err = f"{method} not implemented for {api}"
            raise NotImplementedError(err)

        # Call the API
        data = methods[api](records, method=method)

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
