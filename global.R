library(shiny)
library(DT)
library(stringr)
library(dplyr)
library(RMySQL)
library(pool)
library(DBI)
library(foreach)
library(zip)
library(shinyjs)
source('locuscompare.R')
source('config/config.R')
library(digest)
library(utils)
library(googledrive)
library(googlesheets)
library(promises)
library(future)
plan(multiprocess)

# Variables:
locuscompare_pool = dbPool(
    RMySQL::MySQL(), 
    dbname = "locuscompare",
    host = aws_host,
    username = aws_username,
    password = aws_password,
    minSize = 4,
    maxSize = Inf,
    idleTimeout = 3600000
)

args = list(
    drv = RMySQL::MySQL(),
    dbname = "locuscompare",
    host = aws_host,
    username = aws_username,
    password = aws_password
)