# UI file for LocusCompare
# Boxiang Liu
# 2018-01-01


get_study_list = function(locuscompare_pool,pattern = 'eQTL|GWAS') {
	table_names = dbListTables(locuscompare_pool)
	idx = which(str_detect(table_names, pattern) & !str_detect(table_names,'_trait'))
	table_names = table_names[idx]
	study_category = str_split_fixed(table_names, '_', 2)[, 1]
	study_list = foreach(i = unique(study_category), .combine = c) %do% {
		tmp = list(table_names[study_category == i])
		names(tmp) = i
		return(tmp)
	}
	return(study_list)
}

study_list = get_study_list(locuscompare_pool, pattern = 'eQTL|GWAS')

get_coloc_gwas_list = function(locuscompare_pool){

	gwas = dbGetQuery(
		conn = locuscompare_pool,
		statement = 'select distinct gwas from eCAVIAR;'
	)

	return(gwas)

}

coloc_gwas_list = get_coloc_gwas_list(locuscompare_pool)

shinyUI(fluidPage(
	useShinyjs(),
	extendShinyjs(script = sprintf('%s/init.js',home_dir)),

	###################
	# Loading message #
	###################

	div(
		id = "loading-content",
		h2("Loading...")
	),
	hidden(
		div(
			id = "app-content",
			navbarPage(
				title = 'LocusCompare',
				id = 'navbarPage',

				########################
				# Colocalization panel #
				########################
				
				tabPanel(
					title = 'Colocalization',
					fluidRow(h3('Select Studies')),

					fluidRow(
						column(
							width = 2,
							tags$i(h3('Select GWAS'))
						),
						column(
							width = 5,
							selectInput(
								inputId = 'coloc_gwas',
								label = 'GWAS study',
								choices = c('Choose' = '', coloc_gwas_list),
								width = '100%'
							)
						),
						column(
							width = 5,
							selectizeInput(
								inputId = 'coloc_trait',
								label = 'Trait',
								choices = c('Choose' = ''),
								width = '100%'
							)
						)
					),

					fluidRow(
						column(
							width = 2,
							tags$i(h3('Select eQTL'))
						),
						column(
							width = 10,
							selectizeInput(
								inputId = 'coloc_eqtl',
								label = 'eQTL',
								choices = c('Choose' = ''),
								width = '100%'
							)
						)
					),

					br(),

					fluidRow(
						column(
							width = 12,
							actionButton(
								inputId = 'plot_coloc',
								label = 'Plot colocalization', 
								width = '100%',
								class = "btn-primary"
							)
						)
					),

					br(),

					fluidRow(
						column(
							width = 12,
							shinycssloaders::withSpinner(
								plotOutput(
									outputId = 'coloc',
									click = 'coloc_plot_click',
									dblclick = 'coloc_plot_dblclick',
									brush = brushOpts(id = 'coloc_plot_brush', direction = 'x'),
									height = '250px'
								)
							)
						)
					),

					fluidRow(
						dataTableOutput(outputId = 'coloc_gene')
					),

					br(),

					hidden(
						div(
							id = "plot_locuscompare_button",
							fluidRow(
								column(
									width = 12,
									actionButton(
										inputId = 'coloc_to_locuscompare',
										label = 'Plot LocusCompare', 
										width = '100%',
										class = "btn-primary"
									)
								)
							)
						)
					),

					fluidRow(
						column(
							12,
							shinycssloaders::withSpinner(plotOutput(outputId = 'blank_plot_coloc', height = '1px'))
						)
					)
				),

				#####################
				# Single loci panel #
				#####################

				tabPanel(
					title = 'Single Locus',
					fluidRow(h3('Select Studies')),
					
					# Select study 1:
					fluidRow(h3('Study 1')),

					fluidRow(
						column(2,
								tags$i(h3('Select Published'))
						),
						column(5,
								selectInput(
									inputId = 'study1',
									label = 'Study Name',
									choices = c('Choose' = '', study_list),
									width = "100%"
								)
						),
						column(5,
								selectizeInput(
									inputId = 'trait1',
									label = 'Trait (e.g. phenotype or gene)',
									choices = c('Choose' = ''),
									width = "100%"
								)
						)
					),
					fluidRow(
						column(2,
								tags$i(h3(tags$b('Or'),'Upload'))
						),
						column(5,
								fileInput(
									inputId = 'file1', 
									label = downloadLink(outputId = 'file1_example', label = 'Example file'),
									width = "100%"
								)
						),
						column(5,
								textInput(
									inputId = 'file1_trait', 
									label = 'Label (x-axis)', 
									width = "100%"
								)
						)
					),

					fluidRow(
						span(textOutput(outputId = 'check_study1'),style='color:gray')
					),

					# Select study 2:
					fluidRow(h3('Study 2')),
					fluidRow(
						column(2,
								tags$i(h3('Select Published'))
						),
						column(5,
								selectInput(
									inputId = 'study2',
									label = 'Study Name',
									choices = c('Choose' = '', study_list),
									width = "100%"
								)
						),
						column(5,
								selectizeInput(
									inputId = 'trait2',
									label = 'Trait (e.g. phenotype or gene)',
									choices = c('Choose' = ''),
									width = "100%"
								)
						)
					),
					fluidRow(
						column(2,
								tags$i(h3(tags$b('Or'),'Upload'))
						),
						column(5,
								fileInput(
									inputId = 'file2', 
									label = downloadLink(outputId = 'file2_example', label = 'Example file'),
									width = "100%"
								)
						),
						column(5,
								textInput(
									inputId = 'file2_trait', 
									label = 'Label (y-axis)', 
									width = "100%"
								)
						)
					),


					fluidRow(
						span(textOutput(outputId = 'check_study2'),style='color:gray')
					),

					fluidRow(
						h3('Select a region (maximum 2Mb flanking window)')
					),

					fluidRow(
						column(3,tags$i(h3('SNP'))),
						column(3, textInput(inputId = 'reference_snp', label = 'Reference SNP', placeholder = 'e.g. rs1698683')),
						column(3, numericInput(inputId = 'snp_window', label = 'Flanking Window (Kb)', value = 100))
					),

					fluidRow(
						column(3,tags$i(h3(tags$b('Or'),'Gene'))),
						column(3, textInput(inputId = 'reference_gene', label = 'Reference Gene', placeholder = 'e.g. PHACTR1')),
						column(3, numericInput(inputId = 'gene_window', label = 'Flanking Window (Kb)', value = 100))
					),

					fluidRow(
						column(3, tags$i(h3(tags$b('Or'),'Coordinate'))),
						column(
							3, 
							selectInput(
								inputId = 'chr', 
								label = 'Chromosome', 
								choices = c(Choose='',as.character(c(1:22)),'X')
							)
						),
						column(3, numericInput(inputId = 'start',label = 'Start', value = NULL, min = 1)),
						column(3, numericInput(inputId = 'end', label = 'End', value = NULL, min = 1))
					),

					fluidRow(
						span(textOutput(outputId = 'check_region'),style = 'color:gray')
					),

					br(),

					fluidRow(
						column(
							12,
							actionButton(
								inputId = 'interactive_to_locuscompare', 
								label = 'Plot LocusCompare', 
								width = '100%',
								class = "btn-primary"
							)
						)
					),

					fluidRow(
						column(
							12,
							shinycssloaders::withSpinner(plotOutput(outputId = 'blank_plot', height = '1px'))
						)
					),

					br(),

					fluidRow(
						column(
							4,
							actionButton(
								inputId = 'faq',
								label = 'Frequently Asked Questions',
								width = '100%',
								onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki/1.-FAQ', '_blank')"
							)
						),
						column(
							4,
							actionButton(
								inputId = 'documentation',
								label = 'Documentation',
								width = '100%',
								onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki', '_blank')"
							)
						),
						column(
							4,
							actionButton(
								inputId = 'bugs',
								label = 'Report Bugs',
								width = '100%',
								onclick = "window.open('https://github.com/boxiangliu/locuscompare/issues', '_blank')"
							)
						)
					)
				),

				##########################
				# Interactive plot panel #
				##########################

				tabPanel(
					title = 'Plots',

					fluidRow(
						column(4,
							actionButton(
								inputId = 'back',
								label = '< Back'
							)
						)
					),

					br(),

					fluidRow(
						column(
							6,
							selectizeInput(
								inputId = 'snp',
								label = 'SNP',
								choices = c('Choose' = ''),
								width = "100%"
							)
						),
						column(
							6,
							selectInput(
								inputId = 'population',
								label = 'Population:',
								choices = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'),
								selected = 'EUR',
								width = '100%'
							)
						)
					),

					fluidRow(
						column(
							width = 6,
							shinycssloaders::withSpinner(
								plotOutput(
									outputId = 'locuscompare', 
									height = '500px', 
									click = 'plot_click'
								)
							)
						),
						column(
							width = 6,
							fluidRow(
								column(12,
									shinycssloaders::withSpinner(
										plotOutput(
											outputId = 'locuszoom1',
											click = 'plot_click',
											dblclick = 'plot_dblclick',
											brush = brushOpts(id = 'plot_brush', direction = 'x'),
											height = '250px'
										)
									)
								)
							),
							fluidRow(
								column(12,
								       shinycssloaders::withSpinner(
										plotOutput(
											outputId = 'locuszoom2',
											click = 'plot_click',
											dblclick = 'plot_dblclick',
											brush = brushOpts(id = 'plot_brush', direction = 'x'),
											height = '250px'
										)
									)
								)
							)
						)
					),

					br(),
					br(),

					fluidRow(
						column(
							4,
							h3('SNP information'),
							verbatimTextOutput('snp_info'),
							h3('Download options'),
							numericInput(
								inputId = 'locuscompare_length',
								label = 'LocusCompare side length (inches)',
								value = 8,
								min = 1,
								step = 1,
								width = '100%'
							),
							numericInput(
								inputId = 'locuszoom_height',
								label = 'LocusZoom height (inches)',
								value = 4,
								min = 1,
								step = 1,
								width = '100%'
							),
							numericInput(
								inputId = 'locuszoom_width',
								label = 'LocusZoom width (inches)',
								value = 8,
								min = 1,
								step = 1,
								width = '100%'
							),
							downloadButton('single_download', 'Download',class = 'butt')
						),
						column(
							8,
							fluidRow(
								column(6,
										h3('LD SNPs')
								),
								column(6,
										numericInput(
											inputId = "r2_threshold",
											label = "LD threshold (r2)",
											value = 0.8,
											min = 0.2,
											max = 1,
											step = 0.1,
											width = '100%'
										)
								)
							),
							fluidRow(
								dataTableOutput('ld_snps')
							)
						)
					)
				),

				##############
				# Batch mode #
				##############

				tabPanel(
					title = 'Batch Query',
					div(
						id = 'batch_query',
						fluidRow(h3('Select Studies')),
						# Select study 1:
						fluidRow(h3('Study 1')),
						fluidRow(
							column(2,
									tags$i(h3('Select Published'))
							),
							column(10,
									selectInput(
										inputId = 'batch_study1',
										label = 'Study Name',
										choices = c('Choose' = '', study_list),
										width = "100%"
									)
							)
						),
						fluidRow(
							column(2,
								   tags$i(h3(tags$b('Or'),'Upload'))
							),
							column(10,
									fileInput(
										inputId = 'batch_file1', 
										label = downloadLink(outputId = 'batch_file1_example', label = 'Example file'),
										width = "100%"
									)
							)
						),
						# Select study 2:
						fluidRow(h3('Study 2')),
						fluidRow(
							column(2,
									tags$i(h3('Published'))
							),
							column(10,
									selectInput(
										inputId = 'batch_study2',
										label = 'Study Name',
										choices = c('Choose' = '', study_list),
										width = "100%"
									)
							)
						),
						fluidRow(
							column(2,
								   tags$i(h3(tags$b('Or'),'Upload'))
							),
							column(10,
									fileInput(
										inputId = 'batch_file2', 
										label = downloadLink(outputId = 'batch_file2_example', label = 'Example file'),
										width = "100%"
									)
							)
						),
						# Reference population:
						fluidRow(h3('LD calculation')),
						fluidRow(
							column(2,tags$i(h3('Reference Population'))),
							column(10,
									selectInput(
										inputId = 'batch_population',
										label = 'Population:',
										choices = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'),
										selected = 'EUR',
										width = '100%'
									)
							)
						),
						# Genomic coordinates:
						fluidRow(h3('Genomic regions (maximum 25 regions)')),
						fluidRow(
							column(2,tags$i(h3('Input'))),
							column(10,
									textAreaInput(
										inputId = 'batch_input',
										label = 'Genomic regions',
										width = '100%',
										resize = 'vertical',
										placeholder = 'e.g.\nIL18R1:100kb\nchr9:5215926-7194036\nrs2284033:100kb'
									)
							)
						),
						fluidRow(
							column(2,tags$i(h3(tags$b('Or'),'Upload'))),
							column(10,
									fileInput(
										inputId = 'batch_region', 
										label = downloadLink(outputId = 'batch_region_example', label = 'Example file'),
										width = "100%"
									)
							)
						),
						fluidRow(h3('Job metadata')),
						fluidRow(
							column(2,tags$i(h3('Job name'))),
							column(10,
								textInput(
									inputId = 'batch_job_name', 
									label = 'Job name', 
									placeholder = 'myjob', 
									width = '100%'
								)
							)
						),
						fluidRow(
							column(2,tags$i(h3('Email'))),
							column(10,
    							textInput(
    								inputId = 'batch_job_email',
    								label = 'Email',
    								placeholder = 'e.g. me@domain.com',
    								width = '100%'
    							)
							)
						),
						fluidRow(
						    column(2,tags$i(h3('LocusCompare'))),
						    column(10,
						        numericInput(
    						        inputId = 'batch_locuscompare_length',
    						        label = 'LocusCompare side length (inches)',
    						        value = 8,
    						        min = 1,
    						        step = 1,
    						        width = '100%'
						        )
						    )
						),
						fluidRow(
						    column(2,tags$i(h3('LocusZoom'))),
						    column(10,
						        numericInput(
    						        inputId = 'batch_locuszoom_height',
    						        label = 'LocusZoom height (inches)',
    						        value = 4,
    						        min = 1,
    						        step = 1,
    						        width = '100%'
						        )
						     )
						),
						fluidRow(
						    column(2,''),
						    column(10,
						        numericInput(
    						        inputId = 'batch_locuszoom_width',
    						        label = 'LocusZoom width (inches)',
    						        value = 8,
    						        min = 1,
    						        step = 1,
    						        width = '100%'
						        )
						    )
						),
						fluidRow(
							column(12,
									actionButton(
										inputId = 'submit', 
										label = 'Submit',
										width = '100%',
										class = "btn-primary"
									)
							)
						),
						fluidRow(
							column(
								12,
								textOutput(outputId = 'batch_error')
							)
						),
						br(),
						fluidRow(
							column(
								4,
								actionButton(
									inputId = 'faq',
									label = 'Frequently Asked Questions',
									width = '100%',
									onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki/Frequently-Asked-Questions', '_blank')"
								)
							),
							column(
								4,
								actionButton(
									inputId = 'documentation',
									label = 'Documentation',
									width = '100%',
									onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki', '_blank')"
								)
							),
							column(
								4,
								actionButton(
									inputId = 'bugs',
									label = 'Report Bugs',
									width = '100%',
									onclick = "window.open('https://github.com/boxiangliu/locuscompare/issues', '_blank')"
								)
							)
						)
					),
					shinyjs::hidden(
						div(id = 'batch_query_success',
							h3('Your query has been submitted succesfully!')
						)
					)
				),

				#################
				# Download page #
				#################

				tabPanel(
					title = 'Download',
					fluidRow(
						column(
							width = 12,
							h3('Downloading GWAS and QTL datasets'),
							p('LocusCompare current hosts 50 (TODO: add automatic update) 
							  GWAS datasets and 44 (TODO: add automatic update) eQTL datasets.
							  We have spent a large amount of time finding and curating these 
							  datasets, and make little sense for anyone else to repeat the work.
							  Therefore, we wrote a convenient bash script for you to download 
							  these data.'),
							p('Link to bash script: ',a('https://github.com/mikegloudemans/gwas-download',href='https://github.com/mikegloudemans/gwas-download'))
							)
						),
					fluidRow(
						column(
							width = 12,
							h3('Acknowledgement'),
							p('We are grateful to the individuals and organizations who have made 
							  their data publically available. Below is a list of studies 
							  we included in LocusCompare. If we missed a reference to your 
							  study, please email Boxiang Liu <bliu2@stanford.edu> or 
							  Mike Gloudemans <mgloud@stanford.edu>.'),
							dataTableOutput(outputId = 'study_info')
							)
						)
				),

				###################
				# Contribute page #
				###################

				tabPanel(
					title = 'Contribute',
					div(id = 'form',
						fluidRow(
							column(
								width = 12,
								h3('Share your study'),
								p('Current GWAS datasets are distributed across various websites, making them difficult to find 
								  and download. LocusCompare provides a platform for sharing GWAS and QTL datasets. If you would like
								  to contribute your data to LocusCompare, please fill out the following form. 
								  Alternatively, you can email', tags$a(href="mailto:bliu2@stanford.edu", "Boxiang Liu"), 'or', 
								  tags$a(href="mailto:mgloud@stanford.edu", "Mike Gloudemans"),
								  'to add your study.',tags$b('(* = mandatory)'))
								)
						),
						fluidRow(
							column(
								width = 6,
								textInput(inputId = 'form_trait', label = 'Trait/Molecular Phenotype*', width = '100%'),
								textInput(inputId = 'form_ethnicity', label = 'Ethnicity', width = '100%'),
								textInput(inputId = 'form_sample_size', label = 'Sample Size', width = '100%'),
								textInput(inputId = 'form_author', label = 'Author/Consortium*', width = '100%'),
								fileInput(inputId = 'form_file', label = 'Upload Association Summary Statistics', width = '100%')
							),
							column(
								width = 6,
								textInput(inputId = 'form_year', label = 'Year*', width = '100%'),
								textInput(inputId = 'form_journal', label = 'Journal*', width = '100%'),
								textInput(inputId = 'form_link', label = 'Link to Publication*', width = '100%'),
								textInput(inputId = 'form_download_link', label = 'Link to Dataset', width = '100%'),
								textInput(inputId = 'form_comments',label = 'Comments', width = '100%')
							)
						),
						fluidRow(
							column(
								width = 12,
								actionButton(inputId = 'form_submit', label = 'Submit', class = "btn-primary", width = '100%')
							)
						)
					),
					shinyjs::hidden(
						div(id = 'thankyou_msg',
							h3('Thank you! Your response has been submitted succesfully.'),
							actionLink('submit_another','Share another dataset')
						)
					)
				)
			)
		)
	)
))
