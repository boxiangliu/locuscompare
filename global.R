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

# Variables:
tmp_dir=tempdir()
Sys.chmod(tmp_dir, mode="0777")

locuscompare_pool = dbPool(
	RMySQL::MySQL(), 
	dbname = "locuscompare",
	host = aws_host,
	username = aws_username,
	password = aws_password
)

