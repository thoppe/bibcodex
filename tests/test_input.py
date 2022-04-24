import bibcodex
import pandas as pd
import json
import pytest

# from hypothesis import given, strategies as st, settings

sample_DOIs = [
    "Not a DOI",
    "Also not a DOI",
    "10.1007/s10742-021-00241-z",
    "10.1001/jama.2018.14854",
    None,
]

sample_PMIDs = [
    None,
    "30193304",
    "29516104",
    29516104,  # This shouldn't be a valid PMID
    "Not a PMID",
]


@pytest.fixture(scope="module")
def sample_bibcodex():
    """
    Builds a fixture of a dataframe
    """
    df = pd.DataFrame()
    df.bibcodex.clear()

    df["doi"] = sample_DOIs
    df["pmid"] = sample_PMIDs

    yield df

    # Teardown of fixture here
    df.bibcodex.clear()


@pytest.fixture(scope="module")
def valid_bibcodex(sample_bibcodex):
    """
    Builds a fixture of a single call to PubMed so repeated calls
    are not needed.
    """
    df = sample_bibcodex

    idx = df.bibcodex.valid_pmid_idx()
    df = df[idx].set_index("pmid")
    info = df.bibcodex.download("pubmed")
    df["title"] = info["title"]

    yield df


###########################################################################


def test_cache_size(valid_bibcodex):
    """
    Using the fixture, test that the cache has the size of two after a call.
    """
    df = valid_bibcodex
    assert len(df.bibcodex.pubmed) == 2


def test_cache_disksize(valid_bibcodex):
    """
    Using the fixture, test that the disk size of the cache is non-zero.
    """
    df = valid_bibcodex
    assert df.bibcodex.pubmed.size() > 0


def test_cache_pull(valid_bibcodex):
    """
    Using the fixture, try to redownload items already present.
    """
    df = valid_bibcodex
    df.bibcodex.download("pubmed")


def test_cache_iteration(valid_bibcodex):
    """
    Using the fixture, iterate through the keys and grab them.
    Check that the result is valid json by dumping to a string,
    and checking the key is in that string.
    """
    df = valid_bibcodex
    for key in df.bibcodex.pubmed:
        val = df.bibcodex.pubmed.get(key)
        js = json.dumps(val)
        assert key in js


def test_cache_iteration(valid_bibcodex):
    """
    Using the fixture, check the keys are the same as the iteration.
    """
    df = valid_bibcodex
    for key in df.bibcodex.pubmed.keys():
        assert key in df.bibcodex.pubmed


def test_empty_cache():
    """
    With an empty cache, make sure that the cache size is zero.
    """
    df = pd.DataFrame()
    df.bibcodex.pubmed.clear()
    assert len(df.bibcodex.pubmed) == 0


def test_validation_DOI(sample_bibcodex):
    """
    Check that DOIs are marked correctly.
    """
    df = sample_bibcodex
    info = df.bibcodex.validate()

    assert info["n_rows"] == 5
    assert info["n_doi"] == 4
    assert info["n_doi_invalid"] == 2
    assert info["n_doi_missing"] == 1


def test_validation_PMID(sample_bibcodex):
    """
    Check that PMIDs are marked correctly.
    """
    df = sample_bibcodex
    info = df.bibcodex.validate()

    assert info["n_rows"] == 5
    assert info["n_pmid"] == 4
    assert info["n_pmid_invalid"] == 2
    assert info["n_pmid_missing"] == 1


def test_invalid_API():
    df = pd.DataFrame()

    with pytest.raises(TypeError):
        df.bibcodex.download("NOT A REAL API")
