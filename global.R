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

home_dir='/srv/persistent/bliu2/locuscompare/'


# Variables:
tmp_dir=tempdir()
Sys.chmod(tmp_dir, mode="0777")

locuscompare_db = dbPool(
	RMySQL::MySQL(), 
	dbname = "locuscompare",
	host = "localhost",
	username = "root",
	password = "admin"
)

reference_db = dbPool(
	RMySQL::MySQL(),
	dbname = 'reference',
	host = "localhost",
	username = "root",
	password = "admin"
)

