###############################
# Batch section from server.R #
###############################
#----------------# 
#   Batch mode   #
#----------------#

output$batch_file1_example = downloadHandler(
	filename = function(){return('PHACTR1_Coronary_Heart_Disease_Nikpay_2015_batch.tsv')},
	content = function(file){file.copy(sprintf('%s/data/example/PHACTR1_Coronary_Heart_Disease_Nikpay_2015_batch.tsv',home_dir),file)},
	contentType = 'text/txt'
	)

output$batch_file2_example = downloadHandler(
	filename = function(){return('PHACTR1_all_tissues.tsv')},
	content = function(file){file.copy(sprintf('%s/data/example/PHACTR1_all_tissues.tsv',home_dir),file)},
	contentType = 'text/txt'
	)

output$batch_region_example = downloadHandler(
	filename = function(){return('batch_region.txt')},
	content = function(file){file.copy(sprintf('%s/data/example/batch_region.txt',home_dir),file)},
	contentType = 'text/txt'
	)

batch_observer = reactiveValues(batch_file1 = NULL, batch_file2 = NULL, batch_region_upload = NULL)

# Check batch study 1 is valid:
observeEvent(
	eventExpr = input$batch_file1,
	handlerExpr = {
		batch_observer$batch_file1 = TRUE
	}
)

batch_study1_total = reactive({
	return(isTruthy(input$batch_study1) + isTruthy(batch_observer$batch_file1))
})

output$check_batch_study1 = renderText({
	if (batch_study1_total() == 0){
		return('Please either select or upload Study 1.')
	} else if (batch_study1_total() == 1){
		return('Study 1 has a valid input.')
	} else if (batch_study1_total() == 2){
		return('Please either select or upload a study, but not both.')
	} else {
		return('Please contact Boxiang Liu at jollier.liu@gmail.com')
	}
})

observeEvent(
	eventExpr = input$batch_file1_reset,
	handlerExpr = {
		reset('batch_file1')
		batch_observer$batch_file1 = FALSE 
	}
)

# Check batch study 2 is valid:
observeEvent(
	eventExpr = input$batch_file2,
	handlerExpr = {
		batch_observer$batch_file2 = TRUE
	}
)

batch_study2_total = reactive({
	return(isTruthy(input$batch_study2) + isTruthy(batch_observer$batch_file2))
})

output$check_batch_study2 = renderText({
	if (batch_study2_total() == 0){
		return('Please either select or upload Study 2.')
	} else if (batch_study2_total() == 1){
		return('Study 2 has a valid input.')
	} else if (batch_study2_total() == 2){
		return('Please either select or upload a study, but not both.')
	} else {
		return('Please contact Boxiang Liu at jollier.liu@gmail.com')
	}
})

observeEvent(
	eventExpr = input$batch_file2_reset,
	handlerExpr = {
		reset('batch_file2')
		batch_observer$batch_file2 = FALSE
	}
)

# Check batch region is valid:
observeEvent(
	eventExpr = input$batch_region_upload,
	handlerExpr = {
		batch_observer$batch_region_upload = TRUE
	}
)


batch_region_total = reactive({
	return(isTruthy(input$batch_region_input) + isTruthy(batch_observer$batch_region_upload))
})

output$check_batch_region = renderText({
	if (batch_region_total() == 0){
		return('Please either select or upload genomic regions.')
	} else if (batch_region_total() == 1){
		return('Genomic regions have a valid input.')
	} else if (batch_region_total() == 2){
		return('Please either select or upload genomic regions, but not both.')
	} else {
		return('Please contact Boxiang Liu at jollier.liu@gmail.com')
	}
})

observeEvent(
	eventExpr = input$batch_region_upload_reset,
	handlerExpr = {
		reset('batch_region_upload')
		batch_observer$batch_region_upload = FALSE
	}
)

output$check_batch_job_name = renderText({
	if (isTruthy(input$batch_job_name)){
		return('Valid job name.')
	} else {
		return('Please enter a job name.')
	}
})

output$check_batch_job_email = renderText({

	if (isTruthy(input$batch_job_email)){
	
		if (str_detect(input$batch_job_email,'.+@.+')){

			return('Valid email.')

		} else {

			return('Invalid email.')
		}

	} else {

		return('Please enter an email')

	}
})

batch_submit_ready = reactive({
	batch_study1_total() == 1 & 
		batch_study2_total() == 1 & 
		batch_region_total() == 1 & 
		isTruthy(input$batch_job_name) & 
		isTruthy(str_detect(input$batch_job_email,'.+@.+'))
})

observe({
	shinyjs::toggleState('batch_submit', batch_submit_ready())
})

observeEvent(input$batch_submit,{

	if (isTruthy(input$batch_region_input)){
		coordinate_list = str_split(trimws(input$batch_region_input),'\n')[[1]]
	} else {
		coordinate_list = unlist(fread(input$batch_region_upload$datapath,header=FALSE,sep='\t'))
	}

	coordinate_list = coordinate_list[1:min(25,length(coordinate_list))] # Subset to first 25 because of gmail size limit.
	
	valid_batch_study1_ = isTruthy(input$batch_study1)
	valid_batch_study2_ = isTruthy(input$batch_study2)
	input_ = reactiveValuesToList(input)
	token_ = session$token
	tar_fn = future({batch_query(tmp_dir,coordinate_list,valid_batch_study1_,valid_batch_study2_,input_,token_)})

	link = tar_fn %...>% 
		googledrive::drive_upload(media = ., path = paste0('LocusCompare/Download/',basename(.))) %...>%
		googledrive::drive_share(role = 'reader', type = 'anyone') %...>%
		googledrive::drive_link()

	subject = sprintf('LocusCompare job %s completed on %s',input$batch_job_name,Sys.time())

	msg = link %...>% sprintf('LocusCompare job %s was completed on %s. Download via this link: %s',input$batch_job_name,Sys.time(),.)
	
	msg %...>% send.mail(from = email_username,
		to = input$batch_job_email,
		subject = subject,
		body = .,
		smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = email_username, passwd = email_password, ssl = TRUE),
		authenticate = TRUE,
		send = TRUE)
})

observeEvent(input$batch_submit,{
	shinyjs::hide('batch_query')
	shinyjs::show('batch_query_success')
})


observeEvent(input$submit_another_query,{
	shinyjs::reset('batch_query')
	shinyjs::show('batch_query')
	shinyjs::hide('batch_query_success')
})

###########################
# Batch section from ui.R #
###########################
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
										width = "100%",
										selected = 'eQTL_Lung_GTEx_v7'
									)
							)
						),

						fluidRow(
							column(
								width = 2,
								tags$i(h3(tags$b('Or'),'Upload'))
							),

							column(
								width = 8,
								fileInput(
									inputId = 'batch_file1', 
									label = downloadLink(outputId = 'batch_file1_example', label = 'Example file'),
									width = "100%"
								)
							),

							column(
								width = 2,
								actionButton(
									inputId = 'batch_file1_reset',
									label = 'Reset',
									width = '100%',
									style = 'margin-top:25px')
							)
						),

						fluidRow(
							column(
								width = 12,
								textOutput(outputId = 'check_batch_study1')
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
										width = "100%",
										selected = 'GWAS_2_Hour_Glucose_Saxena_2010'
									)
							)
						),

						fluidRow(
							column(
								width = 2,
								tags$i(h3(tags$b('Or'),'Upload'))
							),

							column(
								width = 8,
								fileInput(
									inputId = 'batch_file2', 
									label = downloadLink(outputId = 'batch_file2_example', label = 'Example file'),
									width = "100%"
								)
							),

							column(
								width = 2,
								actionButton(
									inputId = 'batch_file2_reset',
									label = 'Reset',
									width = '100%',
									style = 'margin-top:25px')
							)
						),

						fluidRow(
							column(
								width = 12,
								textOutput(outputId = 'check_batch_study2')
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
							column(
								width = 4,
								textAreaInput(
									inputId = 'batch_region_input',
									label = 'Genomic regions',
									width = '100%',
									resize = 'vertical',
									placeholder = 'e.g.\nIL18R1:1000000\nchr9:5215926-7194036\nrs2284033:1000000',
									value = 'chr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\nchr1:1000000-2000000\nchr1:2000000-3000000\nchr1:3000000-4000000\n'
								)
							),
							column(
								width = 6,
								helpText(
									"Valid input formats:",
									br(),
									"1. chr<i>:<start>-<end>",
									br(),
									"2. <gene symbol>:<window>",
									br(),
									"3. <rsid>:<window>"
								)
							)
						),

						fluidRow(
							column(2,tags$i(h3(tags$b('Or'),'Upload'))),
							column(
								width = 8,
								fileInput(
									inputId = 'batch_region_upload', 
									label = downloadLink(outputId = 'batch_region_example', label = 'Example file'),
									width = "100%"
								)
							),
							column(
								width = 2,
								actionButton(
									inputId = 'batch_region_upload_reset',
									label = 'Reset',
									width = '100%',
									style = 'margin-top:25px')
							)
						),

						fluidRow(
							column(
								width = 12,
								textOutput(outputId = 'check_batch_region')
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
									width = '100%',
									value = 'test'
								)
							)
						),

						fluidRow(
							column(
								width = 12,
								textOutput(outputId = 'check_batch_job_name')
							)
						),

						fluidRow(
							column(2,tags$i(h3('Email'))),
							column(10,
								textInput(
									inputId = 'batch_job_email',
									label = 'Email',
									placeholder = 'e.g. me@domain.com',
									width = '100%',
									value = 'jollier.liu@gmail.com'
								)
							)
						),

						fluidRow(
							column(
								width = 12,
								textOutput(outputId = 'check_batch_job_email')
							)
						),

						fluidRow(
							column(12,
									actionButton(
										inputId = 'batch_submit', 
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
							h3('Your query has been submitted successfully!'),
							actionLink(
								inputId = 'submit_another_query',
								label = 'Submit another query.',
								width = '100%'
							)
						)
					)
				),
