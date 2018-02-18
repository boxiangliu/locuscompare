# UI file for LocusCompare
# Boxiang Liu
# 2018-01-01


get_study_list = function(locuscompare_db) {
	table_names = dbListTables(locuscompare_db)
	idx = which(!str_detect(table_names, '_traits'))
	table_names = table_names[idx]
	study_category = str_split_fixed(table_names, '_', 2)[, 1]
	study_list = foreach(i = unique(study_category), .combine = c) %do% {
		tmp = list(table_names[study_category == i])
		names(tmp) = i
		return(tmp)
	}
	return(study_list)
}

study_list=get_study_list(locuscompare_db)

shinyUI(fluidPage(
	useShinyjs(),
	extendShinyjs(script = '/srv/persistent/bliu2/locuscompare/init.js'),
	navbarPage(
		title = 'LocusCompare',
		id = 'navbarPage',
		# tabPanel('Instructions',
				 # fluidRow(
					#  column(
					# 	 12,
					# 	 helpText(
					# 		 'LocusCompare is an interactive visualization tool for comparing two genome-wide association datasets. For instance, it can be used to visualize the colocalization between a GWAS dataset to an eQTL dataset, which is particularly useful for fine-mapping studies. The "Visualization" tab has an example dataset to illustrate its usage.'
					# 	 ),
					# 	 hr(),
					# 	 helpText(
					# 		 'On the topright corner is the locuscompare plot. Each axis represent the -log10(P-value) from a dataset. Double-clicking on a point will highlight the selected SNP, and will update the color of other SNPs according to their r2 with the selected SNP. The LD information is calculated using 1000 Genomes phase 3 version 5a. One can update the population-specific LD using the dropdown menu to the left.'
					# 	 ),
					# 	 hr(),
					# 	 helpText(
					# 		 'The two plots below are manhattan plots showing the -log10(P-value) of two datasets. Notice that the highlighted SNPs are synchronized between locuscompare and manhattan plots. Therefore, double-clicking on any plot will update all three.'
					# 	 ),
					# 	 hr(),
					# 	 helpText(
					# 		 'On the bottom is a table that shows selected SNP as well as SNPs passing a given LD threshold. Single-clicking a SNP on any plot will update the table. LD threshold can be set using the "r2 threshold" input.'
					# 	 )
					#  )
				 # )),
		
		# Interactive input panel:
		tabPanel(
			'Interactive Plot',
			fluidRow(h3('Select Studies')),

			# Select study 1:
			fluidRow(h3('Study 1')),
			fluidRow(
				column(2,
					tags$i(h3('Published'))
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
						label = 'Trait',
						choices = c('Choose' = ''),
						width = "100%"
					)
				)
			),
			fluidRow(
				column(2,
					tags$i(h3('Upload'))
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
						label = 'Trait name', 
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
						label = 'Trait',
						choices = c('Choose' = ''),
						width = "100%"
					)
				)
			),
			fluidRow(
				column(2,
					tags$i(h3('Upload'))
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
						label = 'Trait name', 
						width = "100%"
					)
				)
			),
			hr(),
			fluidRow(
				h3('Select a region (choose one)')
				),
			fluidRow(
				column(
					3,tags$i(h3('SNP'))
				),
				column(
					3, textInput(inputId = 'reference_snp', label = 'Reference SNP', placeholder = 'e.g. rs1698683')
				),
				column(
					3, numericInput(inputId = 'snp_window', label = 'Flanking Window (Kb)', value = 100)
				)
			),
			fluidRow(
				column(
					3,tags$i(h3('Gene'))
				),
				column(
					3, textInput(inputId = 'reference_gene', label = 'Reference Gene', placeholder = 'e.g. BRCA1')
				),
				column(
					3, numericInput(inputId = 'gene_window', label = 'Flanking Window (Kb)', value = 100)
				)
			),
			fluidRow(
				column(
					3, tags$i(h3('Coordinate'))
				),
				column(
					3, 
					selectInput(
						inputId = 'chr', 
						label = 'Chromosome', 
						choices = c(Choose='',as.character(c(1:22)),'X')
					)
				),
				column(
					3, numericInput(inputId = 'start',label = 'Start', value = NULL, min = 1)
				),
				column(
					3, numericInput(inputId = 'end', label = 'End', value = NULL, min = 1)
				)
			),
			br(),
			fluidRow(
				column(
					12, 
					actionButton(
						inputId = 'visualize', 
						label = 'Plot!', 
						width = '100%',
						style = 'font-weight: bold'
					)
				)
			),
			fluidRow(
				column(
					12,
					textOutput(outputId = 'interactive_error')
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
						onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki', '_blank')"
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
		
		tabPanel(
			'Plots',
			fluidRow(
				column(4,
					actionButton(
						inputId = 'back',
						label = '< Back',
					)
				)
			),
			fluidRow(
				column(
					4,
					# h3('Instructions'),
					# helpText(
					# 	'1. The table in the bottom shows selected variant and its LD proxies. Single click to select a variant to display in table. Change "r2 threshold" to control number of LD proxies.'
					# ),
					# helpText(
					# 	'2. One variant is highlighted in purple, and other variants are colored according to r2 value with the highlighted variant. Double click to highlight a variant. Select a population upon which the r2 are calculated.'
					# ),
					h3('LD calculation'),
					selectInput(
						'population',
						'Population:',
						choices = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'),
						selected = 'EUR',
						width = '100%'
					),
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
					plotOutput('locuscompare', height = 'auto', click = 'plot_click')
				)
			),
			fluidRow(
				column(
					12,
					plotOutput(
						outputId = 'locuszoom1',
						click = 'plot_click',
						dblclick = 'plot_dblclick',
						brush = brushOpts(id = 'plot_brush', direction = 'x'),
						height = '200px'
					)
				)
			),
			fluidRow(
				column(
					12,
					plotOutput(
						outputId = 'locuszoom2',
						click = 'plot_click',
						dblclick = 'plot_dblclick',
						brush = brushOpts(id = 'plot_brush', direction = 'x'),
						height = '200px'
					)
				)
			),
			fluidRow(h3('LD SNPs')),
			fluidRow(
				column(
					6,
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
				div(
					style = 'overflow-x: scroll', 
					DT::dataTableOutput('ld_snps')
				)
			)
		),

		# Batch mode
		tabPanel(
			'Batch Plot',
			fluidRow(h3('Select Studies')),

			# Select study 1:
			fluidRow(h3('Study 1')),
			fluidRow(
				column(2,
					tags$i(h3('Published'))
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
					tags$i(h3('Upload'))
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
					tags$i(h3('Upload'))
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
			fluidRow(h3('Genomic regions')),
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
				column(2,tags$i(h3('Upload'))),
				column(10,
					fileInput(
						inputId = 'batch_region', 
						label = downloadLink(outputId = 'batch_region_example', label = 'Example file'),
						width = "100%"
					)
				)
			),
			fluidRow(
				column(6,
					actionButton(
						inputId = 'submit', 
						label = 'Submit',
						width = '100%',
						style = 'font-weight:bold'
					)
				),
				tags$head(tags$style(".butt{width:100%;font-weight:bold;}")), # background color and font color
				column(6,
					downloadButton(
						outputId = 'batch_download', 
						label = 'Download',
						class = 'butt'
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
						onclick = "window.open('https://github.com/boxiangliu/locuscompare/wiki', '_blank')"
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
		)
	)
	
))
