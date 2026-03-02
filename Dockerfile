# =============================================================================
# structureR – Docker build
#
# Conforms to the two-stage convention expected by the docker-cache workflow
# at https://github.com/GWMcElfresh/dockerDependencies
#
#   deps    – system libs + all R package dependencies
#   runtime – application code layered on top of deps
# =============================================================================

# ----------------------------------------------------------------------------
# deps stage – install all R and system dependencies
# Optional ARG lets the docker-cache workflow supply a pre-built base image
# for faster incremental builds.
# ----------------------------------------------------------------------------
ARG BASE_IMAGE=rocker/r-ver:4.4.0
FROM ${BASE_IMAGE} AS deps

ARG SKIP_BASE_DEPS=false

# System libraries required for R package compilation.
# Includes dependencies for devtools, Rcpp, and common bioinformatics packages:
#   libgit2-dev        – gert (used by devtools/usethis)
#   libharfbuzz-dev    – textshaping / ragg (used by devtools)
#   libfribidi-dev     – textshaping / ragg
#   libfreetype6-dev   – systemfonts / ragg (ft2build.h)
#   libfontconfig1-dev – systemfonts
#   libpng-dev         – ragg / png
#   libtiff5-dev       – ragg / tiff
#   libjpeg-dev        – ragg / jpeg
#   libwebp-dev        – ragg / webp
#   libcairo2-dev      – Cairo / grDevices
#   libxt-dev          – X11 / grDevices (needed by some R packages)
#   pandoc             – rmarkdown / devtools (vignette building)
#   libglpk-dev        – igraph (Suggests)
#   libhdf5-dev        – hdf5r (Suggests via Seurat)
#   libbz2-dev         – base R / common
#   liblzma-dev        – base R / common
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libfontconfig1-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libwebp-dev \
        libcairo2-dev \
        libxt-dev \
        zlib1g-dev \
        pandoc \
        libglpk-dev \
        libhdf5-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy only the dependency declaration files (keeps cache layers tight)
COPY DESCRIPTION NAMESPACE ./

# Step 1: install remotes, devtools, and testthat explicitly so they are always
# present in the deps image regardless of which DESCRIPTION fields reference them
RUN Rscript -e "install.packages(c('remotes', 'devtools', 'testthat'))"

# Step 2: install all hard dependencies (Imports/Depends) declared in DESCRIPTION
RUN Rscript -e "remotes::install_deps('.', dependencies = 'strong')"

# Step 3: install Suggests. Fail if any fail to ensure robustness for tests.
RUN Rscript -e "options(warn = 2); remotes::install_deps('.', dependencies = TRUE)"

# ----------------------------------------------------------------------------
# runtime stage – copy application code and install the package itself
# ----------------------------------------------------------------------------
FROM deps AS runtime

WORKDIR /workspace
COPY . .

# Build and install the package (compiles Rcpp)
RUN Rscript -e "devtools::install('.', upgrade = 'never')"

CMD ["Rscript", "-e", "library(structureR); cat('structureR loaded successfully\\n')"]
