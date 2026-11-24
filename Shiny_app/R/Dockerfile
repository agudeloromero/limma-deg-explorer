# Limma DEG Explorer — R Shiny
# Build and run from Shiny_app/R (renv.lock must exist; see 00_setup.R):
#   docker build -t limma-deg-explorer-r .
#   docker run --rm -p 3838:3838 limma-deg-explorer-r
#   open http://localhost:3838
# The release workflow sets R_VERSION to the version recorded in renv.lock.

ARG R_VERSION=4.4.3
FROM rocker/shiny:${R_VERSION}

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
        libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /srv/shiny-server/*

WORKDIR /srv/shiny-server

# Restore the exact package versions first (cached unless renv.lock changes),
# using Linux binaries from Posit Package Manager where available
COPY renv.lock /tmp/renv.lock
RUN . /etc/os-release \
 && export RENV_CONFIG_REPOS_OVERRIDE="https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/latest" \
 && R -q -e "install.packages('renv')" \
 && R -q -e "renv::restore(lockfile = '/tmp/renv.lock', library = .libPaths()[1], prompt = FALSE)"

# The app is served at the root URL
COPY app.R helpers.R ./
COPY templates/ templates/
COPY data/ data/

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
