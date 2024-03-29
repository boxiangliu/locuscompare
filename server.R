# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01

options(shiny.maxRequestSize=5*1024^2) 
logs = reactiveValues(user_count = 0, conn_count = 0)
pool_info = dbGetInfo(locuscompare_pool)
message('###############')
message('# APP STARTED #')
message('###############')


#############
# Functions #
#############
select_snp = function(click,merged){
	merged_subset = nearPoints(merged,click,maxpoints=1,threshold = 10)
	snp = merged_subset %>% dplyr::select(rsid) %>% unlist()
	return(snp)
}

select_row = function(click,x){
	subset = nearPoints(x, click, maxpoints = 1)
	return(subset)
}

get_trait=function(study, conn = locuscompare_pool){
	if (study==''){
		trait = ''
	} else if (str_detect(study,'^GWAS_')){
		trait = dbGetQuery(
			conn = conn,
			statement = sprintf(
				"select display_trait 
				from %s_trait;",study)
			)
	} else if (str_detect(study,'QTL')){
		trait = dbGetQuery(
			conn = conn,
			statement = sprintf(
				"select gene_name 
				from %s_trait;",study)
			)
	} else {
		trait = ''
	}
	trait = unname(unlist(trait))
	return(trait)
}

get_study = function(selected_published,study,trait,datapath,coordinate,genome){
	conn = do.call(dbConnect, args)
	on.exit(dbDisconnect(conn))

	if (str_detect(study,'^eQTL_')){
		if (str_detect(trait,'ENSG')){
			res = trait
		} else {
			res = dbGetQuery(
				conn = conn,
				statement = sprintf(
					"select gene_id 
					from gencode_v19_gtex_v6p
					where gene_name = '%s'",
					trait
					)
				)
			trait = res$gene_id[1]
		}
	}  else if (str_detect(study,'QTL'))
        {
                res = dbGetQuery(
                        conn = conn,
                        statement = sprintf(
                                "select gene_id
                                from %s_trait
                                where gene_name = '%s'",
                                study,
                                trait
                                )
                        )
                trait = res$gene_id[1]
        } else if (str_detect(study,'^GWAS')){
                res = dbGetQuery(
                        conn = conn,
                        statement = sprintf(
                                "select trait
                                from %s_trait
                                where display_trait = '%s'",
                                study,
                                trait
                        )
                )
                trait = res$trait[1]
	}
	
	if (selected_published){
		res=dbGetQuery(
			conn = conn,
			statement = sprintf(
				"select t1.rsid, t1.pval 
				from %s as t1 
				join %s as t2
				on t1.rsid = t2.rsid 
				where t1.trait = '%s' 
				and t2.chr = '%s' 
				and t2.pos >= %s 
				and t2.pos <= %s;",
				study,
				genome,
				trait,
				coordinate$chr,
				coordinate$start,
				coordinate$end
				)
			)
	} else {
		res=fread(datapath,header=TRUE,colClasses=c(rsid='character',pval='numeric'))
		shiny::validate(need(all(c('rsid','pval')%in%colnames(res)),'Input file must have columns rsid, pval!'))

		rsid_list=dbGetQuery(
			conn = conn,
			statement = sprintf(
				"select rsid 
				from %s 
				where chr = '%s' 
				and pos >= %s 
				and pos <= %s;",
				genome,
				coordinate$chr,
				coordinate$start,
				coordinate$end
				)
			)
		res = res %>% dplyr::filter(rsid %in% rsid_list$rsid)
	}
	setDT(res)
	return(res)
}

preview_eCAVIAR = function(gwas, trait, conn = locuscompare_pool){

	statement = sprintf(
		"select eqtl, gene_name, clpp
		from eCAVIAR
		where gwas = '%s'
		and trait = '%s'",
		gwas,
		trait
	)
	eCAVIAR_preview = dbGetQuery(
		conn = conn,
		statement = statement
		)
	setDT(eCAVIAR_preview)
	setorder(eCAVIAR_preview, -clpp)
	eCAVIAR_preview_2 = eCAVIAR_preview[,list(Tested = .N, Genes = paste(gene_name, collapse = ',')), by = 'eqtl']
	setnames(eCAVIAR_preview_2,'eqtl','eQTL')
	return(eCAVIAR_preview_2)
}

get_eCAVIAR = function(gwas, trait, eqtl, conn = locuscompare_pool){

	statement = sprintf(
		"select * 
		from eCAVIAR
		where gwas = '%s'
		and trait = '%s'
		and eqtl = '%s'",
		gwas,
		trait,
		eqtl
	)
	eCAVIAR = dbGetQuery(
		conn = conn,
		statement = statement
		)
	setDT(eCAVIAR)
	setnames(eCAVIAR,'chr','chrom')
	if (!str_detect(eCAVIAR$chrom[1],'chr')){
		eCAVIAR[,chrom := paste0('chr',chrom)]
	}
	return(eCAVIAR)
}

epochTime = function() {
	as.integer(Sys.time())
}

humanTime = function() {
	format(Sys.time(), "%Y%m%d-%H%M%OS")
}

saveData = function(data,dir,name) {
	fileName = sprintf("%s_%s_info.csv",
		name,
		humanTime())

	write.csv(x = data, file = file.path(dir, fileName),
		row.names = FALSE, quote = FALSE)
}

fwrite_return = function(x,file,sep){
	fwrite(x=x,file=file,sep=sep)
	return(file)
}

ggsave_return = function(filename,plot,width,height){
	ggsave(filename=filename,plot=plot,width=width,height=height)
	return(filename)
}

get_coloc_trait = function(gwas, conn = locuscompare_pool){

	statement = sprintf("select distinct trait from eCAVIAR where gwas = '%s';", gwas)
	
	result = dbGetQuery(
		conn = conn,
		statement = statement
	)

	return(result$trait)

}

get_coloc_eqtl = function(gwas, trait, conn = locuscompare_pool){

	statement = sprintf("select distinct eqtl from eCAVIAR where gwas = '%s' and trait = '%s';", gwas, trait)
	
	result = dbGetQuery(
		conn = conn,
		statement = statement
	)

	return(result$eqtl)
	
}

shinyServer(function(input, output, session) {


	message('SESSION STARTED.')
	# observe({
	# 	input$preview_coloc
	# 	input$plot_coloc
	# 	input$coloc_to_locuscompare
	# 	input$interactive_to_locuscompare
	# 	message('Number of free connections:', pool_info$numberFreeObjects)
	# 	message('Number of taken connections:', pool_info$numberTakenObjects)
	# })

    onSessionStart = isolate({
        logs$user_count = logs$user_count + 1
        message('Number of users: ',logs$user_count)
    })
    
    onSessionEnded(function(){
        isolate({
            logs$user_count = logs$user_count - 1
            message('SESSION ENDED.')
            message('Number of users: ',logs$user_count)
        })
    })

	#----------------------------#
	# Session-specific variables #
	#----------------------------#

	tmp_dir = paste0(tempdir(),'/',session$token,'/')
	dir.create(tmp_dir,recursive = TRUE)
	Sys.chmod(tmp_dir, mode="0777")

	observeEvent(input$go_to_website,{
		hide(id = "select-genome", anim = TRUE, animType = "fade")    
		show(id = "app-content")
	})

	genome = reactive({
		if (input$genome_assembly == 'GRCh37/hg19') {
			return('tkg_p3v5a_hg19')
		} else {
			return('tkg_p3v5a_hg38')
		}
	})

	transcriptome = reactive({
		if (input$genome_assembly == 'GRCh37/hg19') {
			return('gencode_v19_gtex_v6p')
		} else {
			return('gencode_v19_gtex_v6p_hg38')
		}
	})
	#--------------------#
	#  Interactive mode  #
	#--------------------#

	observeEvent(input$study1,{
		trait1=get_trait(input$study1)
		updateSelectizeInput(session, "trait1", choices = trait1, server = TRUE)
	})

	observeEvent(input$study2,{
		trait2=get_trait(input$study2)
		updateSelectizeInput(session, "trait2", choices = trait2, server = TRUE)
	})

	selected_published_1 = reactive({isTruthy(input$study1) & isTruthy(input$trait1)})
	selected_published_2 = reactive({isTruthy(input$study2) & isTruthy(input$trait2)})

	study1_total = reactive({
		uploaded_new = isTruthy(input$file1) & isTruthy(input$file1_trait)
		return(selected_published_1() + uploaded_new)
	})

	output$check_study1 = renderText({
		if (study1_total() == 0){
			return('Please either select or upload Study 1.')
		} else if (study1_total() == 1){
			return('Study 1 has a valid input.')
		} else if (study1_total() == 2){
			return('Please either select or upload a study, but not both.')
		} else {
			return('Please contact Boxiang Liu at jollier.liu@gmail.com')
		}
	})

	study2_total = reactive({
		uploaded_new = isTruthy(input$file2) & isTruthy(input$file2_trait)
		return(selected_published_2() + uploaded_new)
	})

	output$check_study2 = renderText({
		if (study2_total() == 0){
			return('Please either select or upload Study 2.')
		} else if (study2_total() == 1){
			return('Study 2 has a valid input.')
		} else if (study2_total() == 2){
			return('Please either select or upload a study, but not both.')
		} else {
			return('Please contact Boxiang Liu at jollier.liu@gmail.com')
		}
	})

	selected_snp_region = reactive({isTruthy(input$reference_snp) & isTruthy(input$snp_window)})
	selected_gene_region = reactive({isTruthy(input$reference_gene) & isTruthy(input$gene_window)})
	selected_coordinate = reactive({isTruthy(input$chr) & isTruthy(input$start) & isTruthy(input$end)})

	region_total = reactive({
		return(selected_snp_region() + selected_gene_region() + selected_coordinate())
	})

	output$check_region = renderText({
		if (region_total() == 0){
			return('Please either type in a SNP, a gene, or a coordinate.')
		} else if (region_total() == 1){
			return('Region has a valid input.')
		} else if (region_total() > 1){
			return('Please select only one of SNP, gene, or coordinate.')
		} else {
			return('Please contact Boxiang Liu at jollier.liu@gmail.com')
		}
	})

	interactive_to_locuscompare_ready = reactive({
		study1_total() == 1 & study2_total() == 1 & region_total() == 1
	})

	coloc_to_locuscompare_ready = reactive({
		isTruthy(input$coloc_gwas) & isTruthy(input$coloc_trait) & isTruthy(input$coloc_eqtl) 
	})

	either_to_locuscompare_ready = reactive({
		message(sprintf('Interactive ready? %s',interactive_to_locuscompare_ready()))
		message(sprintf('Colocalization ready? %s',coloc_to_locuscompare_ready()))
		interactive_to_locuscompare_ready() || coloc_to_locuscompare_ready()
	})

	observe({
		shinyjs::toggleState('interactive_to_locuscompare', interactive_to_locuscompare_ready())
	})

	counter = reactiveVal(list(count = 0, from = ''))

	observeEvent(
		eventExpr = {input$interactive_to_locuscompare},
		handlerExpr = {
			new_count = counter()[['count']] + 1
			new_from = 'interactive_to_locuscompare'
			counter(list(count = new_count, from = new_from))
		}
	)

	observeEvent(
		eventExpr = {input$coloc_to_locuscompare},
		handlerExpr = {
			new_count = counter()[['count']] + 1
			new_from = 'coloc_to_locuscompare'
			counter(list(count = new_count, from = new_from))
			message(sprintf('Click count: %s',new_count))
			message(sprintf('Click from: %s',new_from))
		}
	)

	observeEvent(counter(), {
		if (interactive_to_locuscompare_ready()){
			message('Interactive -> LocusCompare')
			showTab(inputId = "navbarPage", target = "Plots", select = TRUE)
		}

		if (coloc_to_locuscompare_ready()){
			shiny::req(coloc_gene_id())
			message('Colocalization -> LocusCompare')
			showTab(inputId = "navbarPage", target = "Plots", select = TRUE)
		} 
	})

	observeEvent(input$back,{

		if (counter()[['from']] == 'interactive_to_locuscompare'){
			message('Back to Interactive')
			showTab(inputId = 'navbarPage', target = 'Single Locus', select = TRUE)
		} else if (counter()[['from']] == 'coloc_to_locuscompare'){
			message('Back to Colocalization')
			showTab(inputId = 'navbarPage', target = 'Colocalization', select = TRUE)
		} 
		hideTab(inputId = "navbarPage", target = "Plots")
		shinyjs::reset('Plots')
	})

	coordinate = eventReactive(counter(),{
		if (counter()[['from']] == 'interactive_to_locuscompare'){
			message('Interactive -> Retrieving coordinates')
			if (selected_snp_region()){

				chr_pos=dbGetQuery(
					conn = locuscompare_pool,
					statement = sprintf('select chr,pos from %s where rsid = "%s" limit 1;',genome(),input$reference_snp)
					)
				shiny::validate(need(nrow(chr_pos)!=0,sprintf('SNP %s not found!',input$reference_snp)))
				shiny::validate(need(nrow(chr_pos)==1,sprintf('SNP %s is not unique!',input$reference_snp)))

				snp_window = min(input$snp_window, 1000) # max window size = 1000kb
				snp_window = max(snp_window, 1) # min window size = 1kb

				res=list(
					chr = chr_pos$chr,
					start = chr_pos$pos - snp_window*500,
					end = chr_pos$pos + snp_window*500
					)

			} else if (selected_gene_region()){

				chr_start_end=dbGetQuery(
					conn = locuscompare_pool,
					statement = sprintf('select chr,start,end from %s where gene_name = "%s" limit 1;', transcriptome(), input$reference_gene)
					)
				shiny::validate(need(nrow(chr_start_end)!=0,sprintf('Gene %s not found!',input$reference_gene)))
				shiny::validate(need(nrow(chr_start_end)==1,sprintf('Gene %s is not unique!',input$reference_gene)))

				gene_window = min(input$gene_window, 1000) # max window size = 1000kb
				gene_window = max(gene_window, 1) # min window size = 1kb

				res=list(
					chr = chr_start_end$chr,
					start = chr_start_end$start - gene_window*500,
					end = chr_start_end$start + gene_window*500
					)

			} else if (selected_coordinate()){

				start_ = max(1, input$start) # min start position is 1

				end_ = min(input$end, start_ + 1e6) # max end position is start + 2e6
				end_ = max(start_, end_) # min end position is start

				res=list(
					chr = input$chr,
					start = start_,
					end = end_
					)

			} else {

				res = NULL

			}

		} else if (counter()[['from']] == 'coloc_to_locuscompare'){
			message('Colocalization -> Retrieving coordinates')
			shiny::validate(need(coloc_gene_id(), label = 'coloc_gene_id'))

			statement = sprintf('select chr,start,end from %s where gene_id = "%s";',transcriptome(), coloc_gene_id())

			chr_start_end = dbGetQuery(
				conn = locuscompare_pool,
				statement = statement
				)

			res = list(
				chr = chr_start_end$chr,
				start = chr_start_end$start - 1e5,
				end = chr_start_end$start + 1e5
				)


		} else {

			shiny::req(FALSE)

		}

		return(res)
	})

	d1 = eventReactive(counter(),{

		shiny::req(either_to_locuscompare_ready())
		message('Getting study 1')

		if (counter()[['from']] == 'interactive_to_locuscompare'){

			selected_published_1_ = selected_published_1()
			input_study1_ = input$study1
			input_trait1_ = input$trait1
			input_file1_datapath_ = input$file1$datapath

		} else if (counter()[['from']] == 'coloc_to_locuscompare'){
			shiny::req(coloc_gene_id())

			selected_published_1_ = TRUE
			input_study1_ = input$coloc_gwas
			input_trait1_ = input$coloc_trait
			input_file1_datapath_ = ''

		} else {

			shiny::req(FALSE)

		}

		get_study(selected_published_1_,input_study1_,input_trait1_,input_file1_datapath_,coordinate(),genome())

	})

	d2 = eventReactive(counter(),{
		
		shiny::req(either_to_locuscompare_ready())
		message('Getting study 2')

		if (counter()[['from']] == 'interactive_to_locuscompare'){
			
			selected_published_2_ = selected_published_2()
			input_study2_ = input$study2
			input_trait2_ = input$trait2
			input_file2_datapath_ = input$file2$datapath

		} else if (counter()[['from']] == 'coloc_to_locuscompare'){

			shiny::req(coloc_gene_id())
			selected_published_2_ = TRUE
			input_study2_ = input$coloc_eqtl
			input_trait2_ = coloc_gene_id()
			input_file2_datapath_ = ''

		} else {

			shiny::req(FALSE)

		}

		get_study(selected_published_2_,input_study2_,input_trait2_,input_file2_datapath_,coordinate(),genome())

	})

	merged = reactive({
		shiny::req(either_to_locuscompare_ready())
		message('Merging studies')

		d1_non_empty = d1() %>% nrow() %>% `>`(0)
		d2_non_empty = d2() %>% nrow() %>% `>`(0)
		shiny::validate(need(d1_non_empty,'No SNP was found in specified region for study 1. Did you specify a correct region?'))
		shiny::validate(need(d1_non_empty,'No SNP was found in specified region for study 2. Did you specify a correct region?'))

		merged = merge(d1(),d2(),by='rsid',suffixes=c('1','2'),all=FALSE)
		merged = merged %>% get_position(.,genome())

		check_overlap = merged %>% nrow() %>% `>`(0)
		shiny::validate(need(check_overlap,'No overlapping SNPs between two studies'))

		merged = merged %>% setDT()
		merged = merged %>% mutate(logp1 = -log10(pval1),logp2 = -log10(pval2))
		message('Finished merging')
		return(merged)
	})

	snp=reactiveVal(value='',label='snp')
	
	observeEvent(merged(),{
		message('Populating SNP field')
		merged() %>% 
			select(rsid) %>% 
			unname () %>%
			unlist() %>% 
			updateSelectizeInput(session, "snp", choices = ., server = TRUE)
	})

	observeEvent(counter(),{
		shiny::req(either_to_locuscompare_ready())
		non_empty_merge = merged() %>% nrow() %>% `>`(0)
		non_empty_merge %>% (
			function(non_empty_merge){
				if (non_empty_merge){
					message('Update selected SNP')
					merged() %>% 
						dplyr::slice(which.min(pval1*pval2)) %>% 
						dplyr::select(rsid) %>% 
						unlist() %>% 
						snp()
				} else {
					snp('')
				}
			}
		)
	})

	observeEvent(input$snp,{
		message('Updating selected SNP')
		snp(input$snp)
	})

	chr=reactive({
		message('Updating chromosome')
		chr=unique(coordinate()$chr)
		shiny::validate(need(length(chr)==1,'Studies must only have one chromosome!'))
		return(chr)
	})

	ld=reactive({
		message('Updating LD')
		retrieve_LD(chr(),snp(),input$population)
	})

	color=reactive({
		message('Updating color')
		merged() %>% 
			dplyr::select(rsid) %>% 
			unlist() %>% unname() %>% 
			assign_color(snp(),ld())
	})

	shape=reactive({
		message('Updating shape')
		merged() %>% assign_shape(snp())
	})

	size=reactive({
		message('Updating size')
		merged() %>% assign_size(snp())
	})

	observeEvent(input$plot_click,{
		selected_snp = merged() %>% select_snp(input$plot_click,.)
		nonempty_snp = selected_snp %>% identical(character(0)) %>% `!`
		if (nonempty_snp){
			selected_snp %>% snp() 
		}

	})

	range=reactiveValues(xmin=NULL,xmax=NULL)

	observeEvent(input$plot_dblclick,{

		brush = input$plot_brush

		if (!is.null(brush)) {
		
			range$xmin = brush$xmin
			range$xmax = brush$xmax

		} else {
		
			range$xmin=NULL
			range$xmax=NULL
		
		}

	})
	
	plot_data=reactive({
		message('Making plot_data')
		if (is.null(range$xmin)){
			plot_data = merged()
		} else {
			plot_data = merged() %>% dplyr::filter(pos<=range$xmax,pos>=range$xmin)
		}
		plot_data = plot_data %>% mutate(label=ifelse(rsid==snp(),rsid,''))
		return(plot_data)
		message('Finished making plot_data')
	})
	
	title1=reactive({

		if (counter()[['from']] == 'interactive_to_locuscompare'){

			if (isTruthy(input$trait1)){
				return(input$trait1)
			} else {
				return(input$file1_trait)
			}

		} else if (counter()[['from']] == 'coloc_to_locuscompare'){

			return(input$coloc_trait)

		} else {

			shiny::req(FALSE)

		}

	})

	title2=reactive({

		if (counter()[['from']] == 'interactive_to_locuscompare'){

			if (isTruthy(input$trait2)){
				return(input$trait2)
			} else {
				return(input$file2_trait)
			}

		} else if (counter()[['from']] == 'coloc_to_locuscompare'){

			return(coloc_gene_id())

		} else {

			shiny::req(FALSE)

		}


	})
	msg = reactive({
		sprintf('%s and %s has no overlapping SNP in %s:%s-%s!',
			title1(), 
			title2(), 
			coordinate()$chr, 
			coordinate()$start, 
			coordinate()$end
	   )
	})

	output$locuscompare = renderPlot({
		nonempty = plot_data() %>% nrow() %>% `>`(0)
		shiny::validate(need(nonempty,msg()))
		message('Making LocusCompare plot')
		p = make_locuscatter(
				merged = plot_data(),
				title1 = title1(),
				title2 = title2(),
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				legend=FALSE
			)
		return(p)
	})
	
	output$locuszoom1 = renderPlot({
		nonempty = plot_data() %>% nrow() %>% `>`(0)
		shiny::validate(need(nonempty,msg()))
		message('Making locuszoom plot 1')
		p = make_locuszoom(
				metal = plot_data(),
				title = title1(),
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				y_string='logp1'
			)
		return(p)
	})
	
	output$locuszoom2 = renderPlot({
		nonempty = plot_data() %>% nrow() %>% `>`(0)
		shiny::validate(need(nonempty,msg()))
		message('Making locuszoom plot 2')
		p = make_locuszoom(
				metal = plot_data(),
				title = title2(),
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				y_string='logp2'
			)
		return(p)
	})
	
	output$blank_plot = renderPlot({
		p = plot_data() %>% {ggplot() + geom_blank()}
		return(p)
	})
	
	output$snp_info = renderText({
		res = dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf('select * from %s where rsid = "%s";',genome(),snp()))

		sprintf(
			'Chromosome: %s\nPosition: %s\nrs ID: %s\nReference SNP: %s\nAlternate SNP: %s\nAllele Frequency: %s\nAFR Frequency: %s\nAMR Frequency: %s\nEAS Frequency: %s\nEUR Frequency: %s\nSAS Frequency: %s',
			res$chr,
			res$pos,
			res$rsid,
			res$ref,
			res$alt,
			res$AF,
			res$AFR_AF,
			res$AMR_AF,
			res$EAS_AF,
			res$EUR_AF,
			res$SAS_AF
			)
	})

	output$ld_snps = renderDataTable({
		shiny::req(coordinate())

		ld_snps=ld() %>%
			dplyr::filter(SNP_A==snp()) %>%
			dplyr::select(rsid=SNP_B,r2=R2)
		
		ld_snps=rbind(data.frame(rsid=snp(),r2=1),ld_snps)
		
		ld_snps_2 = merged() %>%
			dplyr::mutate(pval1_disp = format.pval(pval1), pval2_disp = format.pval(pval2)) %>%
			dplyr::select(rsid,chr,pos,pval1_disp,pval2_disp) %>%
			merge(ld_snps,by='rsid') %>%
			dplyr::rename(rsID = rsid, Chromosome = chr, Position = pos, `P-value (x-axis)` = pval1_disp, `P-value (y-axis)` = pval2_disp, `LD (r2)` = r2) %>%
		    datatable(., selection='none')

		return(ld_snps_2)
	})
	
	output$file1_example = downloadHandler(
		filename = function(){return('PHACTR1_Artery_Coronary.tsv')},
		content = function(file){file.copy(sprintf('%s/data/example/PHACTR1_Artery_Coronary.tsv',home_dir),file)},
		contentType = 'text/txt')

	output$file2_example = downloadHandler(
		filename = function(){return('PHACTR1_Coronary_Heart_Disease_Nikpay_2015.tsv')},
		content = function(file){file.copy(sprintf('%s/data/example/PHACTR1_Coronary_Heart_Disease_Nikpay_2015.tsv',home_dir),file)},
		contentType = 'text/txt')

	output$single_download = downloadHandler(
		filename=function(){return('results.zip')},
		content=function(file){
			data_fn = merged() %>% fwrite_return(paste0(tmp_dir,'/data.tsv'),sep='\t')
			fwrite(ld(),paste0(tmp_dir,'/ld.tsv'),sep='\t')

			locuscompare = make_locuscatter(
					merged = plot_data(),
					title1 = title1(),
					title2 = title2(),
					ld = ld(),
					color = color(),
					shape = shape(),
					size = size(),
					legend=FALSE)
			
			locuszoom1 = make_locuszoom(
					metal = plot_data(),
					title = title1(),
					ld = ld(),
					color = color(),
					shape = shape(),
					size = size(),
					y_string='logp1')
				
			locuszoom2 = make_locuszoom(
					metal = plot_data(),
					title = title2(),
					ld = ld(),
					color = color(),
					shape = shape(),
					size = size(),
					y_string='logp2')

			
			length_ = input$locuscompare_length
			width_ = input$locuszoom_width
			height_ = input$locuszoom_height
			locuscompare_fn = locuscompare %>% ggsave_return(paste0(tmp_dir,'/locuscompare.pdf'),.,width=length_,height=length_)
			locuszoom1_fn = locuszoom1 %>% ggsave_return(paste0(tmp_dir,'/locuszoom1.pdf'),.,width=width_,height=height_)
			locuszoom2_fn = locuszoom2 %>% ggsave_return(paste0(tmp_dir,'/locuszoom2.pdf'),.,width=width_,height=height_)
			c(data_fn,paste0(tmp_dir,'/ld.tsv'),locuscompare_fn,locuszoom1_fn,locuszoom2_fn) %>% utils::zip(file,.,flags = '-j')

		}
		
	)

	#----------------#
	# Colocalization #
	#----------------#

	observeEvent(
		eventExpr = input$coloc_gwas,
		handlerExpr = {
			coloc_trait = get_coloc_trait(input$coloc_gwas)
			updateSelectizeInput(session, "coloc_trait", choices = coloc_trait, server = FALSE)
		}
	)

	observeEvent(
		eventExpr = input$coloc_trait,
		handlerExpr = {
			coloc_eqtl = get_coloc_eqtl(input$coloc_gwas, input$coloc_trait)
			updateSelectizeInput(session, "coloc_eqtl", choices = coloc_eqtl, server = FALSE)
		}
	)

	eCAVIAR_preview = eventReactive(
		eventExpr = input$preview_coloc,
		valueExpr = {
			preview_eCAVIAR(input$coloc_gwas, input$coloc_trait)
		}
	)

	output$coloc_table = renderDataTable(
		expr = eCAVIAR_preview() %>% datatable(.,rownames = FALSE, selection='none',options = list(scrollX = TRUE)),
	)

	eCAVIAR = eventReactive(
		eventExpr = input$plot_coloc,
		valueExpr = {
			get_eCAVIAR(input$coloc_gwas, input$coloc_trait, input$coloc_eqtl)
		}
	)
	


	output$coloc_text = renderText({
		sprintf('[SUMMARY] %s loci passed the threshold (GWAS lead SNP p-value < 5e-8 and eQTL lead SNP p-value < 1e-6).',nrow(eCAVIAR()))
	})
	
	output$coloc_helper_text = renderText({
	    eCAVIAR()
	    '[INSTRUCTIONS]
	    1. Click on a point to see more;
	    2. Drag and double-click to zoom in;
	    3. Double-click again to zoom out.'
	})

	build = 'hg19'
	color1 = 'black'
	color2 = 'black'
	chrom_lengths = manhattan::get_chrom_lengths(build)
	xmax = manhattan::get_total_length(chrom_lengths)
	x_breaks = manhattan::get_x_breaks(chrom_lengths)
	color_map = c(color1,color2)
	names(color_map) = c(color1,color2)

	coloc_range=reactiveValues(x=NULL)

	observeEvent(input$coloc_plot_dblclick, {

		brush = input$coloc_plot_brush

		if (!is.null(brush)) {
		
			coloc_range$x = c(brush$xmin, brush$xmax)

		} else {
		
			coloc_range$x = NULL
		
		}
	
	})

	coloc_plot_data = reactive({

		eCAVIAR() %>% 
			rename(y = clpp) %>%
			manhattan::add_cumulative_pos(.,build) %>%
			manhattan::add_color(., color1 = color1, color2 = color2)

	})

	output$coloc = renderPlot({

		p = coloc_plot_data() %>% 
			{ggplot2::ggplot(.,aes(x=cumulative_pos,y=y,color=color))+
				geom_point(size=5)+
				theme_classic()+
				theme(axis.title = element_text(size=20),axis.text = element_text(size=15)) +
				scale_x_continuous(expand=c(0.01,0),limits=c(0,xmax),breaks=x_breaks,labels=names(x_breaks),name='Chromosome')+
				scale_y_continuous(name='Coloc. Prob.', limits = c(0,NA))+
				scale_color_manual(values=color_map,guide='none')+
				coord_cartesian(xlim = coloc_range$x)}
		return(p)

	})

	selected_row = eventReactive(
		eventExpr = input$coloc_plot_click,
		valueExpr = {
			coloc_plot_data() %>% select_row(input$coloc_plot_click,.)
		}
	)

	coloc_gene_id = reactive({
		selected_row() %>% select(gene_id) %>% unlist() %>% unname()
	})

	output$coloc_gene = renderDataTable({
		selected_row() %>% 
			select(`Gene ID` = gene_id,
				`Gene Symbol` = gene_name, 
				Chromosome = chrom, 
				TSS = pos,
				`GWAS -log10(P)` = logp_gwas,
				`eQTL -log10(P)` = logp_eqtl,
				`Coloc. Prob.` = y) %>% 
			datatable(.,rownames = FALSE, options = list(dom = 't'))
	})

	observe({
		selected_row() %>% (
			function(x){
				if (nrow(x) > 0) {
					shinyjs::show(id = "plot_locuscompare_button")
				} else {
					shinyjs::hide(id = "plot_locuscompare_button")
				}
			}
		)
	})

	observeEvent(
		eventExpr = input$plot_coloc,
		handlerExpr = {
			shinyjs::hide(id = "plot_locuscompare_button")
		}
	)

	output$blank_plot_coloc = renderPlot({
		p = plot_data() %>% {ggplot() + geom_blank()}
		return(p)
	})

	#---------------# 
	# Download page #
	#---------------#

	#sheet_key = googlesheets::gs_key(x='1gq46xlOk674Li50cpv9riZYG7bsfSeZB5qSefa82bR8',lookup=FALSE, verbose = FALSE)
	#list_of_studies = googlesheets::gs_read(sheet_key, verbose = FALSE, col_types = readr::cols())
	#output$study_info = renderDataTable({
	#	DT::datatable(list_of_studies)
	#})

	#-------#
	# Share #
	#-------#

	mandatory_fields = c('form_trait','form_author','form_year','form_journal','form_link')
	observe({
		mandatory_filled = vapply(
			X = mandatory_fields,
			FUN = function(x) isTruthy(input[[x]]),
			FUN.VALUE = logical(1))
		mandatory_filled = all(mandatory_filled)
		shinyjs::toggleState(id = 'form_submit', condition = mandatory_filled)
	})

	description_fields = c('form_trait','form_ethnicity','form_sample_size',
		'form_author','form_year','form_journal','form_link','form_download_link','form_comments')

	formData = eventReactive(input$form_submit,{
		data = sapply(description_fields, function(x) input[[x]])
		data = c(data, file_name = ifelse(isTruthy(input$form_file),input$form_file$name,''))
		data = data.frame(t(data))
		return(data)
	})

	observeEvent(input$form_submit,{
		shinyjs::hide('form')
		shinyjs::show('thankyou_msg')
		file_name = sprintf('%s_%s_%s_%s',input$form_trait,input$form_author,input$form_year,input$form_journal)
		saveData(formData(),contrib_dir,file_name)
		if (isTruthy(input$form_file)) {
			file_path = paste0(contrib_dir,'/',file_name,'_',humanTime(),'_',input$form_file$name)
			file.rename(input$form_file$datapath,file_path)}
		send.mail(from = email_username,
				  to = bosh_email,
				  subject = 'New study uploaded to LocusCompare', 
				  body = sprintf('Path: %s/%s',contrib_dir,file_name), 
				  smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = email_username, passwd = email_password, ssl = TRUE),
				  authenticate = TRUE, 
				  send = TRUE)
	})

})
