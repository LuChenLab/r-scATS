# Installation
## Version Selection
Please choose the appropriate version of scATS based on your R environment:

- `v0.5.6` (Recommended): For users with R >= 4.4.
- `v0.5.5` (Legacy): For users with R < 4.4.


## Requirements
- System Dependencies

Before installing scATS, please ensure that `samtools` is installed and available in your system's PATH.
```bash
conda install -c bioconda samtools
```

- R Dependencies

The following R packages are required for all versions of scATS. You can install them using the following command:

```r
install.packages(c(
    'R.utils', 
    'TTR', 
    'VGAM',
    'prettyunits',
    'httr2',
    'mclust',
    'mixtools'))
```

> **Note:** 💡 For v0.5.5, you must manually install the smoother package from the CRAN archive before installing scATS:
```r
install.packages('https://cran.r-project.org/src/contrib/Archive/smoother/smoother_1.3.tar.gz', repos = NULL, type = 'source')
```

## Installation
To install `scATS`, you have two options: either install directly from GitHub or use the compressed source file.

- Via GitHub (only for `v0.5.6`):
```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("LuChenLab/r-scATS", force = TRUE)

# or
if (!require("remotes")) install.packages("remotes")
remotes::install_github("LuChenLab/r-scATS")
```

- Via Source File (for `v0.5.6` and `v0.5.5`):

Alternatively, you can install `scATS` using the source file downloaded from the [repository](https://github.com/LuChenLab/r-scATS) :
```bash
# Install scATS from a downloaded source file
R CMD INSTALL scATS_0.5.6.tar.gz
```

```r
# or
install.packages("scATS_0.5.6.tar.gz", repos=NULL)

```

## Installation Logs
To ensure reproducibility, we provide comprehensive installation logs for three different installation methods across multiple R versions:
- [Installation Logs for R 4.5.x](https://github.com/LuChenLab/r-scATS/tree/main/logs/R4.5)
- [Installation Logs for R 4.0.x](https://github.com/LuChenLab/r-scATS/tree/main/logs/R4.0)


## Docker Guide

We provide a pre-configured Docker container to eliminate system-specific dependency conflicts.

- Step 1: Install Docker

If you do not have Docker installed, please follow the [official Docker installation guide](https://docs.docker.com/engine/install/).

- Step 2: Run the scATS Container

Pull the pre-built image from the GitHub Container Registry and start an interactive session:

```bash
# Pull the image
docker pull ghcr.io/kayla-xu/scats-evi:r451

# Run the container
docker run -it ghcr.io/kayla-xu/scats-evi:r451
```

- Step 3: Build Docker Locally (Optional)

If you prefer to build the image locally, use the provided Dockerfiles. Building [logs](https://github.com/LuChenLab/r-scATS/tree/main/logs/docker) are also provided for your reference.

```bash
# Build the base environment
docker build --no-cache \
  -t scATS/base:r451 \
  -f Dockerfile_base_r451 . \
  2>&1 | tee Dockerfile_base_r451.log

# Build the final scATS image
docker build --no-cache \
  -t scATS/evi:r451 \
  -f Dockerfile_scats_r451 . \
  2>&1 | tee Dockerfile_scats_r451.log
```