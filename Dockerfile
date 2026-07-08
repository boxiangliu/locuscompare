# LocusCompare single-locus Shiny app.
# Base: rocker/shiny (ships shiny-server + a pinned Posit Package Manager CRAN
# snapshot for reproducible installs).
FROM rocker/shiny:4.4.1

# System libraries. The prebuilt RMySQL binary links against MySQL's client
# library, so we need libmysqlclient21 (provides libmysqlclient.so.21) at runtime;
# libmariadb-dev covers any source builds. The rest are for the usual
# curl/xml/ssl-dependent R packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libmysqlclient21 \
        libmariadb-dev \
        libssl-dev \
        libxml2-dev \
        libcurl4-openssl-dev \
        zlib1g-dev \
        zip \
        git \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages (from the base image's dated snapshot) + the latest locuscomparer.
RUN install2.r --error --skipinstalled \
        shiny shinyjs shinythemes shinycssloaders DT dplyr stringr data.table \
        DBI RMySQL pool ggplot2 cowplot ggrepel markdown remotes \
    && R -e "remotes::install_github('boxiangliu/locuscomparer', upgrade = 'never')" \
    && rm -rf /tmp/downloaded_packages /tmp/*.rds

# App files (see .dockerignore for what is excluded).
RUN rm -rf /srv/shiny-server/*
COPY --chown=shiny:shiny . /srv/shiny-server/locuscompare/
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

ENV APP_HOME=/srv/shiny-server/locuscompare
EXPOSE 3838
# rocker/shiny's default CMD launches shiny-server using the config above.
