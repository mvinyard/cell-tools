
# setup.py

__module_name__ = "setup.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import os
import re
from setuptools import setup
import sys


setup(
    name="cell_tools",
    version="0.0.1",
    python_requires=">3.7.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url="https://github.com/mvinyard/cell-tools",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="cell_tools - Analysis Tools for Single-Cell Data",
    packages=[
        "cell_tools",
        "cell_tools._ATAC",
        "cell_tools._RNA",
        "cell_tools._datasets",
        "cell_tools._multiome",
        "cell_tools._readwrite",
        "cell_tools._utilities",
    ],
    install_requires=[
        "anndata>=0.7.8",
        "cython>=0.29.24",
	"episcanpy>=0.3.2",
	"licorice>=0.0.2",
	"nb-black>=1.0.7",
        "pyranges>=0.0.113",
	"pytrellis>=0.0.1",
	"tqdm>=4.62.3",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
