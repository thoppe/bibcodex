import setuptools
import os

package_name = "bibcodex"
__local__ = os.path.abspath(os.path.dirname(__file__))

f_version = os.path.join(__local__, package_name, "_version.py")
exec(open(f_version).read())

# Get the long description from the relevant file
description = "Access, analyze, and display bibliographic information"

# URL shown on pypi
homepage_url = "https://github.com/thoppe/bibcodex"

# Author information
author_name = "Travis Hoppe"
author_email = "travis.hoppe+{package_name}@gmail.com"

license_name = "CC-SA"

with open("README.md") as FIN:
    long_description = FIN.read()

setuptools.setup(
    name=package_name,
    packages=setuptools.find_packages(),
    # Include package data...
    include_package_data=True,
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=__version__,  # noqa: F821
    # The project's main homepage.
    url=homepage_url,
    # Author details
    author=author_name,
    author_email=author_email,
    # Choose your license
    license=license_name,
    install_requires=[],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5  -Production/Stable
        "Development Status :: 4 - Beta",
        # Indicate who your project is intended for
        "Intended Audience :: Developers",
        "Intended Audience :: Other Audience",
        "Intended Audience :: Science/Research",
        "Topic :: Database",
        "Topic :: Office/Business",
        "Topic :: Scientific/Engineering",
        
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 3",
    ],
    # What does your project relate to?
    keywords=["bibliographic", "publications", "pubmed", "NLP"],
    test_suite="pytest",
    tests_require=["pytest", "hypothesis"],
)
