import diskcache as dc
import brotlicffi
import requests
import pickle
import logging
import pandas as pd
from tqdm import tqdm
from typing import List, Union, Dict


class CachedDownloader:
    """
    Virtual class to cache downloads. Subclass for each data source.
    Data are cached to disk and compressed with Brotli.

    Can only get/save dictionaries or raw text.
    """

    # These must be set by the derived class
    name = None
    datatype = None
    chunksize = None
    api_key = ""

    def __init__(self):
        if self.name is None:
            raise KeyError("Must set name in subclass")

        self.logger = logging.getLogger(f"cache_{self.name}")
        self.cache = dc.Cache(f"cache/{self.name}")
        self.sess = requests.session()

    def __iter__(self):
        yield from self.cache

    def __len__(self):
        # Return the number of items in the cache
        return len(list(self.cache))

    def size(self):
        # Returns the disk size of the cache
        return self.cache.volume()

    def clear(self):
        # Clears the cache
        return self.cache.clear()

    def keys(self):
        return list(self)

    def validate_pmid(self, pmid: int):
        # Raises a TypeError if pmid is not a non-negative int

        err = (
            f"{self.name} expected an non-negative integer for a PMID, "
            f"called with {pmid}"
        )

        if not isinstance(pmid, str):
            raise TypeError("Expect PMIDs to be string data type")

        if isinstance(pmid, str):
            try:
                pmid = int(pmid)
            except (ValueError, TypeError):
                raise TypeError(err)

        if (not pd.api.types.is_number(pmid)) or pmid <= 0 or pd.isnull(pmid):
            raise TypeError(err)

        return True

    def check_pmid(self, pmid: int) -> bool:
        # Returns True/False if the PMID is valid, does not raise an exception
        try:
            self.validate_pmid(pmid)
        except (TypeError, ValueError):
            return False
        return True

    def validate_doi(self, doi: str):
        # Raises a TypeError if doi doesn't match standard

        if not isinstance(doi, str):
            raise TypeError("Expect DOIs to be string data type")

        err = (
            f"{self.name} expected an DOI to start with 10."
            f"called with {doi}"
        )

        if len(doi) < 3 or not doi.startswith("10."):
            raise TypeError(err)

    def check_doi(self, doi: str) -> bool:
        # Returns True/False if the DOI is valid, does not raise an exception
        try:
            self.validate_doi(doi)
        except (TypeError, ValueError):
            return False
        return True

    def download(self, url, params=None, headers=None):
        r = self.sess.get(url, params=params, headers=headers)
        self.logger.info(f"Downloading {r.url}")

        if not r.ok:
            self.logger.error(
                f"Request failed with status code {r.status_code}"
            )

        return r

    def get(self, key):
        key = str(key)

        if key in self.cache:
            self.logger.info(f"CACHE HIT: {key[:100]}")
            val = self.cache[key]
            val = brotlicffi.decompress(val)
            val = pickle.loads(val)
            return val

        return None

    def set(self, key, val):
        key = str(key)

        self.logger.info(f"CACHE SET: {key[:100]}")

        if not isinstance(val, self.datatype) and val is not None:
            raise KeyError(
                f"{self.name} can only cache {self.datatype} objects"
            )

        val = pickle.dumps(val)
        val = brotlicffi.compress(val)
        self.cache[key] = val

    def __contains__(self, key):
        return str(key) in self.cache

    def get_from_PMIDs(self):
        raise NotImplementedError(f"Must define this func for {self.name}")

    def get_from_DOIs(self):
        raise NotImplementedError(f"Must define this func for {self.name}")

    def __call__(
        self, records: Union[int, List], method="pmid"
    ) -> Dict[str, datatype]:
        """
        Downloads (or pulls from cache) a list of pmids or DOIs.
        """
        if method == "pmid":
            if records is not None:
                # Validate the input datatypes and cast to int (if say floats)
                records = [int(p) for p in records if self.validate_pmid(p)]

            result = self.get_from_PMIDs(records)

        elif method == "doi":
            result = self.get_from_DOIs(records)
        else:
            raise NotImplementedError("method must be one of [pmid, doi]")
        return result

    def _chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in tqdm(range(0, len(lst), n), disable=len(lst) == 0):
            yield lst[i : i + n]

    def cached(func):
        """
        If force==True, cache downloader is skipped
        """

        def protect(
            self, keys, force=False, raise_on_missing=False, *args, **kwargs
        ):

            # If a key is passed, turn it into a list
            if not pd.api.types.is_list_like(keys):
                keys = list([keys])

            # Grab the values already known, if forced all items are unknown
            if not force:
                vals = [self.get(k) for k in keys]
            else:
                vals = [None for k in keys]

            # Find a list of the missing items and compute them
            mkeys = [k for (k, v) in zip(keys, vals) if v is None]

            # Call the function in appropriate sized chunks
            mvals = {}
            for chunk in self._chunks(mkeys, self.chunksize):
                result = func(self, chunk, *args, **kwargs)
                mvals.update(result)

            # Add any found values to the cache
            for k in mkeys:

                if str(k) not in mvals:
                    logging.warning(f"{k} not found in {self.name}")
                    mvals[str(k)] = self.datatype()

                self.set(k, mvals[str(k)])

            vals = [
                mvals[str(k)] if v is None else v for k, v in zip(keys, vals)
            ]

            return vals

        return protect
