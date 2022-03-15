import diskcache as dc
import brotlicffi
import requests
import pickle
from wasabi import msg


class CachedDownloader:
    """
    Virtual class to cache downloads. Subclass for each data source.
    Data are cached to disk and compressed with Brotli.

    Can only get/save dictionaries or raw text.
    """

    name = None
    datatype = dict
    has_multi_cache = False

    def __init__(self):
        if self.name is None:
            raise KeyError("Must set name in subclass")

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

        if not isinstance(pmid, int) or pmid <= 0:
            err = f"{self.name} expected an non-negative integer (PMID), called with {pmid}"
            raise TypeError(err)

    def download(self, url, params=None):
        r = self.sess.get(url, params=params)

        msg.info(f"Downloading {r.url}")

        if not r.ok:
            msg.fail(f"Request failed with status code {r.status_code}")

        return r

    def get(self, key):
        key = str(key)

        if key in self.cache:
            val = self.cache[key]
            val = brotlicffi.decompress(val)
            val = pickle.loads(val)
            return val

        return None

    def set(self, key, val):
        key = str(key)

        if not isinstance(val, self.datatype):
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

        def protect(self, key, *args, **kwargs):

            is_force = "force" in kwargs and kwargs.pop("force") == True

            if not isinstance(key, tuple):

                val = self.get(key)

                if val is not None and not is_force:
                    return val

                val = func(self, key, *args, **kwargs)
                self.set(key, val)

                return val

            else:
                if not self.has_multi_cache:
                    err = f"Multiple caching not implemented for {self.name}"
                    raise NotImplementedError(err)

                # Get the items that exist already
                vals = [self.get(x) for x in key]

                # Find a list of the missing items
                if is_force:
                    missing_keys = key
                else:
                    missing_keys = [k for (k, v) in zip(key, vals) if v is None]
                missing_vals = func(self, missing_keys, is_multi=True, *args, **kwargs)

                # Add any found values to the cache
                for k, v in missing_vals.items():
                    if v is not None:
                        self.set(k, v)

                vals = [missing_vals[k] if v is None else v for k, v in zip(key, vals)]

                return vals

        return protect


class RemoteItemNotFound(Exception):
    """Raised when a requested item is not found on the remote server"""

    def __init__(self, key, database_name):
        msg = f"{key} not found in {database_name}"
        super().__init__(msg)
