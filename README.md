# ![cell-tools_logo](docs/imgs/cell-tools.logo.svg)

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/cell-tools.svg)](https://pypi.python.org/pypi/cell-tools/)
[![PyPI version](https://badge.fury.io/py/cell-tools.svg)](https://badge.fury.io/py/cell-tools)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**cell-tools** is a **Single-Cell Data Analysis Toolkit** with the following areas of focus and development:

1. scRNA-seq data
2. scATAC-seq data
3. Multiomic (mostly scRNA-seq + scATAC-seq) data

Funtion development will likely focus on scATAC-seq given that there are so many python-implemented tools developed already for scRNA-seq. However, wrappers that enable batch-mode QC and pipeline-oriented solutions have been designed and will eventually be added.

### Installation
<!-- **Install with `pip`**:
```BASH
pip install cell-tools
``` -->

**Install the development package**:
```BASH
# (1) clone this repository
git clone https://github.com/mvinyard/cell-tools.git

# (2) install the local project in editable mode
cd ./perturb-tools; pip install -e .
```
