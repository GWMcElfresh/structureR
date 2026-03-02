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

# System libraries required for R package compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy only the dependency declaration files (keeps cache layers tight)
COPY DESCRIPTION NAMESPACE ./

# Install R package dependencies declared in DESCRIPTION
# devtools is added here so the test command can call devtools::test()
RUN Rscript -e " \
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('remotes', 'devtools'), quiet = TRUE); \
    remotes::install_deps('.', dependencies = TRUE, quiet = TRUE) \
    "

# ----------------------------------------------------------------------------
# runtime stage – copy application code and install the package itself
# ----------------------------------------------------------------------------
FROM deps AS runtime

WORKDIR /workspace
COPY . .

# Build and install the package (compiles Rcpp)
RUN Rscript -e " \
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    devtools::install('.', upgrade = 'never', quiet = FALSE) \
    "

CMD ["Rscript", "-e", "library(structureR); cat('structureR loaded successfully\\n')"]
