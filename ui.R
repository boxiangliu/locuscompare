# UI file for LocusCompare
# Boxiang Liu
# 2018-01-01


get_study_list = function(locuscompare_pool) {
	table_names = dbListTables(locuscompare_pool)
	idx = which(str_detect(table_names, 'eQTL|GWAS') & !str_detect(table_names,'_trait'))
	table_names = table_names[idx]
	study_category = str_split_fixed(table_names, '_', 2)[, 1]
	study_list = foreach(i = unique(study_category), .combine = c) %do% {
		tmp = list(table_names[study_category == i])
		names(tmp) = i
		return(tmp)
	}
	return(study_list)
}

study_list=get_study_list(locuscompare_pool)

shinyUI(fluidPage(
	useShinyjs(),
	extendShinyjs(script = sprintf('%s/init.js',home_dir)),
	navbarPage(
		title = 'LocusCompare',
		id = 'navbarPage',
		# Batch mode
		tabPanel(
			'Batch Plot',
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
					tags$i(h3('Or Upload'))
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
					tags$i(h3('OR Upload'))
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
				column(2,tags$i(h3('Or Upload'))),
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
						style = 'font-weight:bold; background-color:#ADD8E6;'
					)
				),
				tags$head(tags$style(".butt{width:100%;font-weight:bold;background-color:#ADD8E6}")), # background color and font color
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
