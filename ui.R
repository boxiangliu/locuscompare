# UI file for LocusCompare
# Boxiang Liu
# 2018-01-01

library(shiny)
library(DT)
library(pool)
library(DBI)
library(stringr)
library(foreach)

locuscompare_db <- dbPool(
    RMySQL::MySQL(), 
    dbname = "locuscompare",
    host = "rds-mysql-locuscompare.cbhpzvkzr3rc.us-west-1.rds.amazonaws.com",
    username = "admin",
    password = "12345678"
)

get_study_list=function(locuscompare_db){
    table_names=dbListTables(locuscompare_db)
    study_category=str_split_fixed(table_names,'_',2)[,1]
    study_list=foreach(i=unique(study_category),.combine=c)%do%{
        tmp=list(table_names[study_category==i])
        names(tmp)=i
        return(tmp)
    }
    return(study_list)
}

study_list=get_study_list(locuscompare_db)

shinyUI(fluidPage(

  # Application title
  navbarPage('LocusCompare',
             tabPanel('Instructions',
                      fluidRow(
                          column(12,
                              helpText('LocusCompare is an interactive visualization tool for comparing two genome-wide association datasets. For instance, it can be used to visualize the colocalization between a GWAS dataset to an eQTL dataset, which is particularly useful for fine-mapping studies. The "Visualization" tab has an example dataset to illustrate its usage.'),
                              hr(),
                              helpText('On the topright corner is the locuscompare plot. Each axis represent the -log10(P-value) from a dataset. Double-clicking on a point will highlight the selected SNP, and will update the color of other SNPs according to their r2 with the selected SNP. The LD information is calculated using 1000 Genomes phase 3 version 5a. One can update the population-specific LD using the dropdown menu to the left.'),
                              hr(),
                              helpText('The two plots below are manhattan plots showing the -log10(P-value) of two datasets. Notice that the highlighted SNPs are synchronized between locuscompare and manhattan plots. Therefore, double-clicking on any plot will update all three.'),
                              hr(),
                              helpText('On the bottom is a table that shows selected SNP as well as SNPs passing a given LD threshold. Single-clicking a SNP on any plot will update the table. LD threshold can be set using the "r2 threshold" input.')))),
             
             tabPanel('Browse',
                      fluidRow(h3('Select Studies')),
                      fluidRow(column(6,selectInput(inputId='study1',label='Study 1',choices=c('Choose'='',study_list)),
                               helpText('or'),
                               fileInput(inputId='upload_study1',label=NULL)),
                               column(6,helpText('Select a study from the dropdown menu or upload your own study. The file should have tab-delimited with four columns and must include the same header. Example'),
                                      tableOutput('exampleInput'))),
                      fluidRow(column(6,selectInput(inputId='study2',label='Study 2',choices=c('Choose'='',study_list)),
                               helpText('or'),
                               fileInput(inputId='upload_study2',label=NULL))),
                      hr(),
                      fluidRow(h3('Select Genomic Region')),
                      fluidRow(column(6,textInput(inputId='locus',label='Locus')),
                               column(6,helpText(br(),'Example:',br(),'PHACTR1 or chr6:11716606-13716606'))),
                      fluidRow(column(6,actionButton(inputId='visualize',label='Plot!'))),
                      hr(),
                      fluidRow(h3('Batch Input')),
                      fluidRow(column(6,textAreaInput(inputId='batch_input',label='Batch input',resize='both',height='200px')),
                               column(6,helpText(br(),'Each line is a gene name or genomic coordinate. Example: ',br(),br(),'PHACTR1',br(),'chr15:66063763-68063763',br(),'SORT1'))),
                      fluidRow(column(6,actionButton(inputId='submit',label='submit')))),
             
             tabPanel('Upload',
                      fluidRow(h3('Upload Studies')),
                      fluidRow(fileInput(inputId='upload_study1',label='Study 1')),
                      fluidRow(fileInput(inputId='upload_study2',label='Study 2'))),
             
             tabPanel('Visualize',
                      fluidRow(column(4,
                                      h3('Instructions'),
                                      helpText('1. The table in the bottom shows selected variant and its LD proxies. Single click to select a variant to display in table. Change "r2 threshold" to control number of LD proxies.'),
                                      helpText('2. One variant is highlighted in purple, and other variants are colored according to r2 value with the highlighted variant. Double click to highlight a variant. Select a population upon which the r2 are calculated.'),
                                      numericInput("r2_threshold", "r2 threshold:", 0.8, min = 0, max = 1),
                                      selectInput('population','Population:',choices=c('AFR','AMR','EAS','EUR','SAS'),selected='EUR'),
                                      downloadButton('download','Download')),
                               column(8,plotOutput('locuscompare',height='auto',click='plot_click'))),
                      fluidRow(column(12,plotOutput('locuszoom1',click='plot_click',height='200px'))),
                      fluidRow(column(12,plotOutput('locuszoom2',click='plot_click',height='200px'))),
                      fluidRow(div(style = 'overflow-x: scroll', DT::dataTableOutput('snp_info')))))

))
