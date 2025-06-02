# Installation

```{note}
This guide explains how to install **LANDMARK**, including dependencies and post-installation checks.
```

## Prerequisites

Before installing **LANDMARK**, ensure you have:

- **Python 3.10 or newer** installed. Check your version with:

  ```bash
  python --version
  ```

- **pip** and (optionally) **virtual environment tools** installed:

  ```bash
  python -m ensurepip --default-pip
  python -m pip install --upgrade pip virtualenv
  ```

It is recommended to install LANDMARK in a fresh virtual environment (with [venv](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/en/latest/)).

## Installation Steps

### Clone the Repository

Install **LANDMARK** by cloning the GitHub repository.

#### HTTP Method

```bash
git clone https://github.com/IGNF/landmark.git
cd landmark
pip install -e .
```

#### SSH Method

```bash
git clone git@github.com:IGNF/landmark.git
cd landmark
pip install -e .
```

## Post-Installation Verification

To confirm the package is correctly installed, run:

```bash
landmark --help
```

You should see the help message with all available command-line options.

---

**Note:** If you encounter any installation issues related to dependencies (such as rasterio, shapely, or geopandas), please ensure your environment uses compatible versions of Python and pip. Using a virtual environment is strongly recommended to avoid dependency conflicts.
