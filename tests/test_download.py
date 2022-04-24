import bibcodex
import pandas as pd
from hypothesis import given, strategies as st, settings
import pytest

sample_DOIs = [
    "10.1001/jama.2017.18444",
    "10.1001/jama.2018.0126",
    "10.1001/jama.2018.0708",
]

sample_PMIDs = [
    "29411024",
    "29372233",
    "29509866",
]

sample_titles = [
    "Progressive Massive Fibrosis in Coal Miners From 3 Clinics in Virginia",
    "Birth Defects Potentially Related to Zika Virus Infection During Pregnancy in the United States",
    "Association Between Estimated Cumulative Vaccine Antigen Exposure Through the First 23 Months of Life and Non–Vaccine-Targeted Infections From 24 Through 47 Months of Age",
]


def fixture_sample_dataframe():
    df = pd.DataFrame()
    df["doi"] = sample_DOIs
    df["pmid"] = sample_PMIDs
    df["title"] = sample_titles
    return df


def test_doi2pmid_with_doi():
    """
    Resolve DOIs from PubMed's ID_converter.
    Check that they match known values.
    """
    df = fixture_sample_dataframe().set_index("doi")
    df.bibcodex.clear()

    info = df.bibcodex.download("doi2pmid")
    assert (info["pmid"] == df["pmid"]).all()


def test_pubmed_with_pmid():
    """
    Resolve publications from PubMed.
    Check if DOIs match known values.
    """
    df = fixture_sample_dataframe().set_index("pmid")
    df.bibcodex.clear()

    info = df.bibcodex.download("pubmed")
    assert (info["doi"] == df["doi"]).all()


def test_icite_with_pmid():
    """
    Resolve publications from iCite.
    Check if DOIs match known values.
    """
    df = fixture_sample_dataframe().set_index("pmid")
    df.bibcodex.clear()

    info = df.bibcodex.download("icite")
    assert (info["doi"] == df["doi"]).all()


def test_semanticScholar_with_pmid():
    """
    Resolve publications from Semantic Scholar using PMIDs.
    Due to API restrictions, only check one value.
    Check if DOIs match known value.
    """
    df = fixture_sample_dataframe().set_index("pmid")[:1]
    df.bibcodex.clear()

    info = df.bibcodex.download("semanticScholar")
    assert (info["doi"] == df["doi"]).all()


def test_semanticScholar_with_doi():
    """
    Resolve publications from Semantic Scholar using DOIs.
    Due to API restrictions, only check one value.

    Semantic Scholar does not return PMIDs
    so check if title matches known value.
    """
    df = fixture_sample_dataframe().set_index("doi")[:1]
    df.bibcodex.clear()

    info = df.bibcodex.download("semanticScholar")
    assert (info["title"] == df["title"]).all()


"""
Old test for reference
@given(st.integers())
def test_NAN_from_int(z):
    #Create NAN from a int. Check that it evaluates to the float rep.
    assert z == NAN(z)
"""
