# UI for the LocusCompare single-locus app.
# Boxiang Liu

shinyUI(fluidPage(
    useShinyjs(),

    ##################
    # Landing screen #
    ##################
    div(
        id = "select-genome",
        h2("LocusCompare"),
        p("Interactively visualize the colocalization of two association datasets ",
          "(e.g. a GWAS and an eQTL) at a single locus."),
        fluidRow(
            column(
                width = 12,
                selectInput(
                    inputId = 'genome_assembly',
                    label = 'Genome assembly',
                    choices = c('GRCh37/hg19', 'GRCh38/hg38'),
                    selected = 'GRCh37/hg19'
                )
            )
        ),
        fluidRow(
            column(
                width = 12,
                actionButton(
                    inputId = 'go_to_website',
                    label = 'Go to website',
                    class = "btn-primary"
                )
            )
        )
    ),

    hidden(
        div(
            id = "app-content",
            navbarPage(
                title = 'LocusCompare',
                id = 'navbarPage',
                theme = shinytheme('cerulean'),

                #############
                # About tab #
                #############
                tabPanel(
                    title = 'About',
                    fluidRow(
                        column(width = 1),
                        column(width = 10, includeMarkdown('about.md')),
                        column(width = 1)
                    )
                ),

                ####################
                # Single Locus tab #
                ####################
                tabPanel(
                    title = 'Single Locus',
                    fluidRow(h3('Upload two association datasets')),
                    helpText('Each file must be tab- or comma-delimited with two columns: ',
                             tags$code('rsid'), ' and ', tags$code('pval'), '.'),

                    # Study 1
                    fluidRow(h3('Study 1')),
                    fluidRow(
                        column(6,
                            fileInput(
                                inputId = 'file1',
                                label = downloadLink(outputId = 'file1_example', label = 'Example file'),
                                placeholder = 'Upload your file (5MB max)',
                                width = "100%"
                            )
                        ),
                        column(6,
                            textInput(
                                inputId = 'file1_trait',
                                label = 'Label (x-axis)',
                                value = '',
                                width = "100%"
                            )
                        )
                    ),
                    fluidRow(span(textOutput(outputId = 'check_study1'), style = 'color:gray')),

                    # Study 2
                    fluidRow(h3('Study 2')),
                    fluidRow(
                        column(6,
                            fileInput(
                                inputId = 'file2',
                                label = downloadLink(outputId = 'file2_example', label = 'Example file'),
                                placeholder = 'Upload your file (5MB max)',
                                width = "100%"
                            )
                        ),
                        column(6,
                            textInput(
                                inputId = 'file2_trait',
                                label = 'Label (y-axis)',
                                value = '',
                                width = "100%"
                            )
                        )
                    ),
                    fluidRow(span(textOutput(outputId = 'check_study2'), style = 'color:gray')),

                    # Region
                    fluidRow(h3('Select a region (maximum window = 1Mb)')),
                    fluidRow(
                        column(3, tags$i(h3('SNP'))),
                        column(3, textInput(inputId = 'reference_snp', label = 'Reference SNP', placeholder = 'e.g. rs9349379')),
                        column(3, numericInput(inputId = 'snp_window', label = 'Flanking Window (Kb)', value = 100, min = 1, max = 1000))
                    ),
                    fluidRow(
                        column(3, tags$i(h3(tags$b('Or'), 'Gene'))),
                        column(3, textInput(inputId = 'reference_gene', label = 'Reference Gene', placeholder = 'e.g. PHACTR1')),
                        column(3, numericInput(inputId = 'gene_window', label = 'Flanking Window (Kb)', value = 100, min = 1, max = 1000))
                    ),
                    fluidRow(
                        column(3, tags$i(h3(tags$b('Or'), 'Coordinate'))),
                        column(3, selectInput(inputId = 'chr', label = 'Chromosome', choices = c(Choose = '', as.character(1:22), 'X'))),
                        column(3, numericInput(inputId = 'start', label = 'Start', value = NULL, min = 1)),
                        column(3, numericInput(inputId = 'end', label = 'End', value = NULL, min = 1))
                    ),
                    fluidRow(span(textOutput(outputId = 'check_region'), style = 'color:gray')),

                    br(),
                    fluidRow(
                        column(12,
                            actionButton(
                                inputId = 'interactive_to_locuscompare',
                                label = 'Plot LocusCompare',
                                width = '100%',
                                class = "btn-primary"
                            )
                        )
                    ),
                    fluidRow(
                        column(12, shinycssloaders::withSpinner(plotOutput(outputId = 'blank_plot', height = '1px')))
                    ),
                    br(),
                    fluidRow(
                        column(4, actionButton(inputId = 'faq', label = 'Frequently Asked Questions', width = '100%',
                            onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki/FAQ', '_blank')")),
                        column(4, actionButton(inputId = 'documentation', label = 'Documentation', width = '100%',
                            onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki', '_blank')")),
                        column(4, actionButton(inputId = 'bugs', label = 'Report Bugs', width = '100%',
                            onclick = "window.open('https://github.com/boxiangliu/locuscompare/issues', '_blank')"))
                    )
                ),

                #############
                # Plots tab #
                #############
                tabPanel(
                    title = 'Plots',
                    fluidRow(column(4, actionButton(inputId = 'back', label = '< Back'))),
                    br(),
                    fluidRow(
                        column(4,
                            selectizeInput(inputId = 'snp', label = 'SNP', choices = c('Choose' = ''), width = "100%")
                        ),
                        column(4,
                            selectInput(inputId = 'population', label = 'Population',
                                choices = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'), selected = 'EUR', width = '100%')
                        )
                    ),
                    fluidRow(
                        column(12,
                            helpText('[INSTRUCTIONS]
                                1. Click on a point to select the corresponding SNP;
                                2. Drag and double-click on the plots on the right to zoom in;
                                3. Double-click again to zoom out.')
                        )
                    ),
                    fluidRow(
                        column(6,
                            shinycssloaders::withSpinner(
                                plotOutput(outputId = 'locuscompare', height = '500px', click = 'plot_click')
                            )
                        ),
                        column(6,
                            fluidRow(column(12,
                                shinycssloaders::withSpinner(
                                    plotOutput(outputId = 'locuszoom1', dblclick = 'plot_dblclick',
                                        brush = brushOpts(id = 'plot_brush', direction = 'x', resetOnNew = TRUE), height = '250px')
                                )
                            )),
                            fluidRow(column(12,
                                shinycssloaders::withSpinner(
                                    plotOutput(outputId = 'locuszoom2', dblclick = 'plot_dblclick',
                                        brush = brushOpts(id = 'plot_brush', direction = 'x', resetOnNew = TRUE), height = '250px')
                                )
                            ))
                        )
                    ),
                    br(), br(),
                    fluidRow(
                        column(4,
                            h3('SNP information'),
                            verbatimTextOutput('snp_info'),
                            h3('Download options'),
                            numericInput(inputId = 'locuscompare_length', label = 'LocusCompare side length (inches)', value = 8, min = 1, step = 1, width = '100%'),
                            numericInput(inputId = 'locuszoom_height', label = 'LocusZoom height (inches)', value = 4, min = 1, step = 1, width = '100%'),
                            numericInput(inputId = 'locuszoom_width', label = 'LocusZoom width (inches)', value = 8, min = 1, step = 1, width = '100%'),
                            downloadButton('single_download', 'Download', class = 'butt')
                        ),
                        column(8,
                            fluidRow(column(12, h3('LD SNPs'))),
                            fluidRow(dataTableOutput('ld_snps'))
                        )
                    )
                )
            )
        )
    )
))
