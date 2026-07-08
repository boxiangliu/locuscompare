# Global setup for the LocusCompare single-locus Shiny app.
#
# Plotting, LD retrieval, and SNP-position lookup are provided by the `locuscomparer`
# package (v1.1.0), which reads the shared 1000 Genomes reference tables from the
# `colotool` database. The app keeps its own connection pool for the reference
# queries it makes directly (region variants, SNP info, gene coordinates).
# Boxiang Liu

library(shiny)
library(shinyjs)
library(shinythemes)
library(shinycssloaders)
library(DT)
library(dplyr)
library(stringr)
library(data.table)
library(DBI)
library(RMySQL)
library(pool)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(locuscomparer)

source('config/config.R')

# Connection pool to colotool for the app's own reference queries.
locuscompare_pool <- dbPool(
    drv = RMySQL::MySQL(),
    dbname = db_name,
    host = db_host,
    port = db_port,
    username = db_user,
    password = db_pass,
    minSize = 1,
    maxSize = 8,
    idleTimeout = 60000
)

onStop(function() {
    poolClose(locuscompare_pool)
})
