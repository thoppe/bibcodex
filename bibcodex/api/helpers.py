import diskcache as dc
import brotlicffi
import requests
import pickle
import logging
import pandas as pd


class CachedDownloader:
    """
    Virtual class to cache downloads. Subclass for each data source.
    Data are cached to disk and compressed with Brotli.

    Can only get/save dictionaries or raw text.
    """

    name = None
    datatype = None

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

        if not pd.api.types.is_integer(pmid) or pmid <= 0:
            err = f"{self.name} expected an non-negative integer for a PMID, called with {pmid}"
            raise TypeError(err)

    def download(self, url, params=None):
        r = self.sess.get(url, params=params)
        self.logger.info(f"Downloading {r.url}")

        if not r.ok:
            self.logger.error(f"Request failed with status code {r.status_code}")

        return r

    def get(self, key):
        key = str(key)

        if key in self.cache:
            self.logger.info(f"CACHE HIT: {key}")
            val = self.cache[key]
            val = brotlicffi.decompress(val)
            val = pickle.loads(val)
            return val

        return None

    def set(self, key, val):
        self.logger.info(f"CACHE SET: {key}")

        key = str(key)

        if not isinstance(val, self.datatype) and val is not None:
            raise KeyError(f"{self.name} can only cache {self.datatype} objects")

        val = pickle.dumps(val)
        val = brotlicffi.compress(val)
        self.cache[key] = val

    def __contains__(self, key):
        return str(key) in self.cache

    def cached(func):
        """
        If force==True, cache downloader is skipped
        """

        def protect(self, keys, force=False, raise_on_missing=False, *args, **kwargs):

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

            # Only call the download function if there are missing keys
            if mkeys:
                mvals = func(self, mkeys, *args, **kwargs)
            else:
                mvals = []

            # Add any found values to the cache
            for k in mkeys:

                if str(k) not in mvals:
                    logging.warning(f"{k} not found in {self.name}")
                    mvals[str(k)] = self.datatype()

                self.set(k, mvals[str(k)])

            vals = [mvals[str(k)] if v is None else v for k, v in zip(keys, vals)]

            return vals

        return protect


class RemoteItemNotFound(Exception):
    """Raised when a requested item is not found on the remote server"""

    def __init__(self, key, database_name):
        msg = f"{key} not found in {database_name}"
        super().__init__(msg)
