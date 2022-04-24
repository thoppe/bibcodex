import bibcodex
import pandas as pd
from hypothesis import given, strategies as st, settings
import pytest

sample_DOIs = [
    "Not a DOI",
    "Also not a DOI",
    "10.1007/s10742-021-00241-z",
    None,
]

sample_PMIDs = [
    None,
    "29372233",
    29372233,  # This shouldn't be a valid PMID
    "Not a PMID",
]


def fixture_sample_dataframe():
    df = pd.DataFrame()
    df["doi"] = sample_DOIs
    df["pmid"] = sample_PMIDs
    return df


def test_validation_DOI():
    """
    Check that DOIs are marked correctly.
    """
    df = fixture_sample_dataframe()
    info = df.bibcodex.validate()

    assert info["n_rows"] == 4
    assert info["n_doi"] == 3
    assert info["n_doi_invalid"] == 2
    assert info["n_doi_missing"] == 1


def test_validation_PMID():
    """
    Check that PMIDs are marked correctly.
    """
    df = fixture_sample_dataframe()
    info = df.bibcodex.validate()

    assert info["n_rows"] == 4
    assert info["n_pmid"] == 3
    assert info["n_pmid_invalid"] == 2
    assert info["n_pmid_missing"] == 1
