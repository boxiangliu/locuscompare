# UI file for LocusCompare
# Boxiang Liu
# 2018-01-01

library(shiny)
library(DT)

shinyUI(fluidPage(

  # Application title
  navbarPage('LocusCompare',
             tabPanel('Instructions',fluidRow(column(12,
                                              helpText('LocusCompare is an interactive visualization tool for comparing two genome-wide association datasets. For instance, it can be used to visualize the colocalization between a GWAS dataset to an eQTL dataset, which is particularly useful for fine-mapping studies. The "Visualization" tab has an example dataset to illustrate its usage.'),
                                              hr(),
                                              helpText('On the topright corner is the locuscompare plot. Each axis represent the -log10(P-value) from a dataset. Double-clicking on a point will highlight the selected SNP, and will update the color of other SNPs according to their r2 with the selected SNP. The LD information is calculated using 1000 Genomes phase 3 version 5a. One can update the population-specific LD using the dropdown menu to the left.'),
                                              hr(),
                                              helpText('The two plots below are manhattan plots showing the -log10(P-value) of two datasets. Notice that the highlighted SNPs are synchronized between locuscompare and manhattan plots. Therefore, double-clicking on any plot will update all three.'),
                                              hr(),
                                              helpText('On the bottom is a table that shows selected SNP as well as SNPs passing a given LD threshold. Single-clicking a SNP on any plot will update the table. LD threshold can be set using the "r2 threshold" input.')))),
             tabPanel('Visualize',
                      fluidRow(column(4,
                                      h3('Instructions'),
                                      helpText('1. The table in the bottom shows selected variant and its LD proxies. Single click to select a variant to display in table. Change "r2 threshold" to control number of LD proxies.'),
                                      helpText('2. One variant is highlighted in purple, and other variants are colored according to r2 value with the highlighted variant. Double click to highlight a variant. Select a population upon which the r2 are calculated.'),
                                      numericInput("r2_threshold", "r2 threshold:", 0.8, min = 0, max = 1),
                                      selectInput('population','Population:',choices=c('AFR','AMR','EAS','EUR','SAS'),selected='EUR')),

                               column(8,plotOutput('locuscompare',height='auto',click='plot_click',dblclick = "plot_dblclick"))),
                      
                      fluidRow(column(12,plotOutput('locuszoom1',click='plot_click',
                                                    dblclick = "plot_dblclick",height='200px'))),
                      fluidRow(column(12,plotOutput('locuszoom2',click='plot_click',
                                                    dblclick = "plot_dblclick",height='200px'))),
                      fluidRow(div(style = 'overflow-x: scroll', DT::dataTableOutput('info')))))

))
