# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01

options(shiny.maxRequestSize=100*1024^2) 

# Functions:
select_snp=function(click,merged){
	merged_subset=nearPoints(merged,click,maxpoints=1)
	snp=merged_subset[,rsid]
	return(snp)
}

convert_unit=function(x){
	x = tolower(x)
	shiny::validate(need(str_detect(x,'[0-9]+(kb|mb)*'),sprintf('%s is not recognized. Try e.g. 100kb.',x)))

	if (str_detect(x,'kb')){
		y = as.integer(str_replace(x,'kb',''))*1e3
	} else if (str_detect(x,'mb')){
		y = as.integer(str_replace(x,'mb',''))*1e6
	} else {
		y= as.integer(x)
	}
	return(y)
}

parse_coordinate=function(coordinate){
	coordinate=trimws(coordinate)
	stopifnot(length(coordinate)==1)

	if (str_detect(coordinate,'^rs|^ss')){ # SNP:window format
		split_coordinate = str_split_fixed(coordinate,':',2)
		reference_snp = split_coordinate[,1]
		snp_window = convert_unit(split_coordinate[,2])

		chr_pos=dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf('select chr,pos from tkg_p3v5a where rsid = "%s";',reference_snp)
			)
		shiny::validate(need(nrow(chr_pos)!=0,sprintf('SNP %s not found!',reference_snp)))
		shiny::validate(need(nrow(chr_pos)==1,sprintf('SNP %s is not unique!',reference_snp)))

		res=list(
			chr = chr_pos$chr,
			start = chr_pos$pos - snp_window,
			end = chr_pos$pos + snp_window
			)
	} else if (str_detect(coordinate,'^chr')){ # chr:start-end format
		split_coordinate = str_split_fixed(coordinate,':',2)

		chr=str_replace(split_coordinate[,1],'chr','')
		pos=split_coordinate[,2]

		split_pos=str_split_fixed(pos,'-',2)
		start=as.integer(split_pos[,1])
		end=as.integer(split_pos[,2])

		res=list(
			chr = chr,
			start = start,
			end = end
			)
	} else {# gene:window format
		split_coordinate = str_split_fixed(coordinate,':',2)
		reference_gene = split_coordinate[,1]
		gene_window = convert_unit(split_coordinate[,2])

		chr_start_end=dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf('select chr,start,end from gencode_v19_gtex_v6p where gene_name = "%s";',reference_gene)
			)

		shiny::validate(need(nrow(chr_start_end)!=0,sprintf('Gene %s not found!',reference_gene)))
		shiny::validate(need(nrow(chr_start_end)==1,sprintf('Gene %s is not unique!',reference_gene)))

		res=list(
			chr = chr_start_end$chr,
			start = chr_start_end$start - gene_window,
			end = chr_start_end$end + gene_window
			)
	}
	return(res)
}

get_trait=function(study, conn = locuscompare_pool){
	if (study==''){
		trait = ''
	} else if (str_detect(study,'^GWAS_')){
		trait = dbGetQuery(
			conn = conn,
			statement = sprintf(
				"select distinct trait 
				from %s;",study)
			)
	} else if (str_detect(study,'^eQTL_')){
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

get_study = function(valid_study,study,trait,datapath,coordinate){
	if (str_detect(study,'^eQTL')){
		res = dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf(
				"select gene_id 
				from gencode_v19_gtex_v6p
				where gene_name = '%s'",
				trait
				)
			)
		trait = res$gene_id[1]
	}

	if (valid_study){
		res=dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf(
				"select t1.rsid, t1.pval 
				from %s as t1 
				join tkg_p3v5a as t2 
				on t1.rsid = t2.rsid 
				where t1.trait = '%s' 
				and t2.chr = '%s' 
				and t2.pos >= %s 
				and t2.pos <= %s;",
				study,
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
			conn = locuscompare_pool,
			statement = sprintf(
				"select rsid 
				from tkg_p3v5a 
				where chr = '%s' 
				and pos >= %s 
				and pos <= %s;",
				coordinate$chr,
				coordinate$start,
				coordinate$end
				)
			)
		res=res[rsid%in%rsid_list$rsid,]
	}
	setDT(res)
	return(res)
}

get_batch_study = function(valid_study,study,datapath,coordinate){
	if (valid_study){
		res=dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf(
				"select t1.trait, t1.rsid, t1.pval 
				from %s as t1 
				join tkg_p3v5a as t2 
				on t1.rsid = t2.rsid  
				and t2.chr = '%s' 
				and t2.pos >= %s 
				and t2.pos <= %s;",
				study,
				coordinate$chr,
				coordinate$start,
				coordinate$end
				)
			)
	} else {
		res=fread(datapath,header=TRUE,colClasses=c(trait='character',rsid='character',chr='character',pos='integer',pval='numeric'))
		shiny::validate(need(all(c('trait','rsid','pval')%in%colnames(res)),'Input file must have columns trait, rsid, pval!'))

		rsid_list=dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf(
				"select rsid 
				from tkg_p3v5a 
				where chr = '%s' 
				and pos >= %s 
				and pos <= %s;",
				coordinate$chr,
				coordinate$start,
				coordinate$end
				)
			)
		res=res[rsid%in%rsid_list$rsid,]
	}
	setDT(res)
	return(res)
}

shinyServer(function(input, output, session) {
	#---------------------#
	#   Interactive mode  #
	#---------------------#
	observeEvent(input$study1,{
		trait1=get_trait(input$study1)
		updateSelectizeInput(session, "trait1", choices = trait1, server = TRUE)
	})
	
	observeEvent(input$study2,{
		trait2=get_trait(input$study2)
		updateSelectizeInput(session, "trait2", choices = trait2, server = TRUE)
	})

	valid_study1 = eventReactive(input$visualize,{isTruthy(input$study1) & isTruthy(input$trait1)})
	valid_study2 = eventReactive(input$visualize,{isTruthy(input$study2) & isTruthy(input$trait2)})

	valid_file1 = eventReactive(input$visualize,{isTruthy(input$file1)})
	valid_file2 = eventReactive(input$visualize,{isTruthy(input$file2)})

	valid_snp_region = eventReactive(input$visualize,{isTruthy(input$reference_snp) & isTruthy(input$snp_window)})
	valid_gene_region = eventReactive(input$visualize,{isTruthy(input$reference_gene) & isTruthy(input$gene_window)})
	valid_coordinate = eventReactive(input$visualize,{isTruthy(input$chr) & isTruthy(input$start) & isTruthy(input$end)})
	
	output$interactive_error = renderText({
		shiny::validate(need(valid_study1() | valid_file1(),'Please provide study 1!'))
		shiny::validate(need(!(valid_study1() & valid_file1()),'Please select or upload study 1, but not both!'))
		
		shiny::validate(need(valid_study2() | valid_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_study2() & valid_file2()),'Please select or upload study 2, but not both!'))

		shiny::validate(need(any(valid_snp_region(),valid_gene_region(),valid_coordinate()),'Please provide a region!'))
		shiny::validate(need({valid_snp_region()+valid_gene_region()+valid_coordinate()==1},'Please only provide one of SNP, gene, or coordinate!'))
	})

	observeEvent(input$visualize, {
		shiny::validate(need(valid_study1() | valid_file1(),'Please provide study 1!'))
		shiny::validate(need(!(valid_study1() & valid_file1()),'Please select or upload study 1, but not both!'))
		
		shiny::validate(need(valid_study2() | valid_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_study2() & valid_file2()),'Please select or upload study 2, but not both!'))

		shiny::validate(need(any(valid_snp_region(),valid_gene_region(),valid_coordinate()),'Please provide a region!'))
		shiny::validate(need({valid_snp_region()+valid_gene_region()+valid_coordinate()==1},'Please only provide one of SNP, gene, or coordinate!'))

		showTab(inputId = "navbarPage", target = "Plots", select = TRUE)
	})

	observeEvent(input$back,{
		hideTab(inputId = "navbarPage", target = "Plots")
	})

	coordinate = eventReactive(input$visualize,{
		shiny::validate(need(any(valid_snp_region(),valid_gene_region(),valid_coordinate()),'Please provide a region!'))
		shiny::validate(need({valid_snp_region()+valid_gene_region()+valid_coordinate()==1},'Please only provide one of SNP, gene, or coordinate!'))

		if (valid_snp_region()){
			chr_pos=dbGetQuery(
				conn = locuscompare_pool,
				statement = sprintf('select chr,pos from tkg_p3v5a where rsid = "%s";',input$reference_snp)
				)
			shiny::validate(need(nrow(chr_pos)!=0,sprintf('SNP %s not found!',input$reference_snp)))
			shiny::validate(need(nrow(chr_pos)==1,sprintf('SNP %s is not unique!',input$reference_snp)))

			res=list(
				chr = chr_pos$chr,
				start = chr_pos$pos - input$snp_window*1e3,
				end = chr_pos$pos + input$snp_window*1e3
				)
		}
		if (valid_gene_region()){
			chr_start_end=dbGetQuery(
				conn = locuscompare_pool,
				statement = sprintf('select chr,start,end from gencode_v19_gtex_v6p where gene_name = "%s";',input$reference_gene)
				)
			shiny::validate(need(nrow(chr_start_end)!=0,sprintf('Gene %s not found!',input$reference_gene)))
			shiny::validate(need(nrow(chr_start_end)==1,sprintf('Gene %s is not unique!',input$reference_gene)))

			res=list(
				chr = chr_start_end$chr,
				start = chr_start_end$start - input$gene_window*1e3,
				end = chr_start_end$end + input$gene_window*1e3
				)
		}
		if (valid_coordinate()){
			res=list(
				chr = input$chr,
				start = input$start,
				end = input$end
				)
		}
		return(res)
	})

    
	d1 = eventReactive(input$visualize,{
		shiny::validate(need(valid_study1() | valid_file1(),'Please provide study 1!'))
		shiny::validate(need(!(valid_study1() & valid_file1()),'Please select or upload study 1, but not both!'))
		get_study(valid_study1(),input$study1,input$trait1,input$file1$datapath,coordinate())
	})

	d2 = eventReactive(input$visualize,{
		shiny::validate(need(valid_study2() | valid_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_study2() & valid_file2()),'Please select or upload study 2, but not both!'))
		get_study(valid_study2(),input$study2,input$trait2,input$file2$datapath,coordinate())
	})

	merged=reactive({
		merged=merge(d1(),d2(),by='rsid',suffixes=c('1','2'),all=FALSE)
		merged=get_position(merged)
		shiny::validate(need(nrow(merged)>0,'No overlapping SNPs between two studies'))
		setDT(merged)
		merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
		return(merged)
	})
	
	snp=reactiveVal(value=NULL,label='snp')
	
	observeEvent(input$visualize,{
		snp(merged()[which.min(pval1*pval2),rsid])
	})

	chr=reactive({
		chr=unique(coordinate()$chr)
		shiny::validate(need(length(chr)==1,'Studies must only have one chromosome!'))
		return(chr)
	})

	ld=reactive({
		retrieve_LD(chr(),snp(),input$population)
	})

	color=reactive({assign_color(merged()$rsid,snp(),ld())})
	shape=reactive({assign_shape(merged(),snp())})
	size=reactive({assign_size(merged(),snp())})
	
	observeEvent(input$plot_click,{
		selected_snp=select_snp(input$plot_click,merged())
		if (!identical(selected_snp,character(0))){
			snp(selected_snp)    
		}
	})

	range=reactiveValues(xmin=NULL,xmax=NULL)
	observeEvent(input$plot_dblclick,{
		brush = input$plot_brush
		if (!is.null(brush)){
			range$xmin=brush$xmin
			range$xmax=brush$xmax
		} else {
			range$xmin=NULL
			range$xmax=NULL
		}
	})
	
	plot_data=reactive({
		if (is.null(range$xmin)){
			plot_data=merged()
		} else {
			plot_data=merged()[pos<=range$xmax & pos>=range$xmin]
		}
		
		plot_data[,label:=ifelse(rsid==snp(),rsid,'')]
		return(plot_data)
	})
	
	title1=reactive({
		if (isTruthy(input$trait1)){
			return(input$trait1)
		} else {
			return(input$file1_trait)
		}
	})

	title2=reactive({
		if (isTruthy(input$trait2)){
			return(input$trait2)
		} else {
			return(input$file2_trait)
		}
	})

	output$locuscompare = renderPlot({
		make_locuscatter(
			merged = plot_data(),
			title1 = title1(),
			title2 = title2(),
			ld = ld(),
			color = color(),
			shape = shape(),
			size = size(),
			legend=FALSE
			)
	},height=function(){
		session$clientData$output_locuscompare_width
	})
	
	output$locuszoom1 = renderPlot({
		make_locuszoom(
			metal = plot_data()[,list(rsid,chr,pos,logp1,label)],
			title = title1(),
			ld = ld(),
			color = color(),
			shape = shape(),
			size = size(),
			y_string='logp1')
	})
	
	output$locuszoom2 = renderPlot({
		make_locuszoom(
			metal = plot_data()[,list(rsid,chr,pos,logp2,label)],
			title = title2(),
			ld = ld(),
			color = color(),
			shape = shape(),
			size = size(),
			y_string='logp2')
	})
	
	output$snp_info = renderText({
		res=dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf('select * from tkg_p3v5a where rsid = "%s";',snp())
			)
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

	output$ld_snps = DT::renderDataTable({
		ld_snps=ld()[SNP_A==snp(),][R2>=input$r2_threshold,list(rsid=SNP_B,r2=R2)]
		ld_snps=rbind(data.table(rsid=snp(),r2=1),ld_snps)
		tmp=merged()[rsid%in%ld_snps$rsid,list(rsid,chr,pos,pval1,pval2)]
		snp_info=merge(tmp,ld_snps,by='rsid')
		DT::datatable(snp_info)
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
			owd=setwd(tmp_dir)
			on.exit({setwd(owd)})
			
			fwrite(merged(),'data.tsv',sep='\t')
			fwrite(ld(),'ld.tsv',sep='\t')

			locuscompare = make_locuscatter(
				merged = plot_data(),
				title1 = title1(),
				title2 = title2(),
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				legend=FALSE
				)

			locuszoom1 = make_locuszoom(
				metal = plot_data()[,list(rsid,chr,pos,logp1,label)],
				title = title1(),
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				y_string='logp1'
				)

			locuszoom2 = make_locuszoom(
				metal = plot_data()[,list(rsid,chr,pos,logp2,label)],
				title = title2(),
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				y_string='logp2'
				)

			ggsave('locuscompare.pdf',locuscompare,width=input$locuscompare_length,height=input$locuscompare_length)
			ggsave('locuszoom1.pdf',locuszoom1,width=input$locuszoom_width,height=input$locuszoom_height)
			ggsave('locuszoom2.pdf',locuszoom2,width=input$locuszoom_width,height=input$locuszoom_height)
			
			zip(file,c('data.tsv','ld.tsv','locuscompare.pdf','locuszoom1.pdf','locuszoom2.pdf'))
		}
		
	)

	#----------------# 
	# 	Batch mode   #
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

	valid_batch_study1 = eventReactive(input$submit,{isTruthy(input$batch_study1)})
	valid_batch_study2 = eventReactive(input$submit,{isTruthy(input$batch_study2)})
	
	valid_batch_file1 = eventReactive(input$submit,{isTruthy(input$batch_file1)})
	valid_batch_file2 = eventReactive(input$submit,{isTruthy(input$batch_file2)})
	
	valid_batch_input = eventReactive(input$submit,{isTruthy(input$batch_input)})
	valid_batch_region = eventReactive(input$submit,{isTruthy(input$batch_region)})

	output$batch_error = renderText({
		shiny::validate(need(valid_batch_study1() | valid_batch_file1(),'Please provide study 1!'))
		shiny::validate(need(!(valid_batch_study1() & valid_batch_file1()),'Please select or upload study 1, but not both!'))
		
		shiny::validate(need(valid_batch_study2() | valid_batch_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_batch_study2() & valid_batch_file2()),'Please select or upload study 2, but not both!'))

		shiny::validate(need(any(valid_batch_input(),valid_batch_region()),'Please provide a list of regions!'))
		shiny::validate(need({valid_batch_input()+valid_batch_region()==1},'Please either input or upload a list of regions, but not both!'))

		return('All inputs are given. Ready to submit!')
	})

	observeEvent(input$submit,{
		shiny::validate(need(valid_batch_study1() | valid_batch_file1(),'Please provide study 1!'))
		shiny::validate(need(!(valid_batch_study1() & valid_batch_file1()),'Please select or upload study 1, but not both!'))
		
		shiny::validate(need(valid_batch_study2() | valid_batch_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_batch_study2() & valid_batch_file2()),'Please select or upload study 2, but not both!'))

		shiny::validate(need(any(valid_batch_input(),valid_batch_region()),'Please provide a list of regions!'))
		shiny::validate(need({valid_batch_input()+valid_batch_region()==1},'Please either input or upload a list of regions, but not both!'))

		owd=setwd(tmp_dir)
		on.exit({setwd(owd)})

		if (valid_batch_input()){
			coordinate_list = str_split(trimws(input$batch_input),'\n')[[1]]
		} else {
			coordinate_list = unlist(fread(input$batch_region$datapath,header=FALSE,sep='\t'))
		}
		
		progress = shiny::Progress$new()
		on.exit(progress$close())
		progress$set(message = "Making plot", value = 0)
		n = length(coordinate_list)

		for (coordinate in coordinate_list){
			if (!dir.exists(paste0(tmp_dir,'/',coordinate))){
				dir.create(paste0(tmp_dir,'/',coordinate),recursive=TRUE)
			}
			
			parsed_coordinate=parse_coordinate(coordinate)


			d1 = get_batch_study(
				valid_study = valid_batch_study1(),
				study = input$batch_study1,
				datapath = input$batch_file1$datapath,
				coordinate = parsed_coordinate
				)

			d2 = get_batch_study(
				valid_study = valid_batch_study2(),
				study = input$batch_study2,
				datapath = input$batch_file2$datapath,
				coordinate = parsed_coordinate
				)

			trait1_list=unique(d1$trait)
			trait2_list=unique(d2$trait)
			n1_trait=length(trait1_list)
			n2_trait=length(trait2_list)

			for (trait1 in trait1_list){
				for (trait2 in trait2_list){
					d1_trait = d1[trait==trait1,list(rsid,pval)]
					d2_trait = d2[trait==trait2,list(rsid,pval)]

					progress$inc(1/n/n1_trait/n2_trait, detail = paste("coordinate: ", coordinate))

					merged=merge(d1_trait,d2_trait,by='rsid',suffixes=c('1','2'),all=FALSE)
					merged=get_position(merged)

					if (nrow(merged)==0) {
						warning(sprintf('%s and %s do not overlap!',trait1, trait2))
						next
					}
					setDT(merged)
					merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
					
					snp=merged[which.min(pval1*pval2),rsid]
					chr=unique(parsed_coordinate$chr)
					if (length(chr)!=1){
						warning(sprintf('%s is not a legal chromosome!',chr))
					}

					ld=retrieve_LD(chr,snp,input$batch_population)

					color=assign_color(merged$rsid,snp,ld)
					shape=assign_shape(merged,snp)
					size=assign_size(merged,snp)
					
					plot_data=copy(merged)
					plot_data[,label:=ifelse(rsid==snp,rsid,'')]
					
					p1=make_locuscatter(
						merged = plot_data,
						title1 = trait1,
						title2 = trait2,
						ld = ld,
						color = color,
						shape = shape,
						size = size,
						legend=FALSE)
					
					p2=make_locuszoom(
						metal = plot_data[,list(rsid,chr,pos,logp1,label)],
						title = trait1,
						ld = ld,
						color = color,
						shape = shape,
						size = size,
						y_string='logp1')
					
					
					p3=make_locuszoom(
						metal = plot_data[,list(rsid,chr,pos,logp2,label)],
						title = trait2,
						ld = ld,
						color = color,
						shape = shape,
						size = size,
						y_string='logp2')
					
					cowplot::save_plot(
						filename = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-locuscompare.pdf'),
						plot = p1)
					
					cowplot::save_plot(
						filename = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-locuszoom1.pdf'),
						plot = p2)
					
					cowplot::save_plot(
						filename = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-locuszoom2.pdf'),
						plot = p3)
					
					data.table::fwrite(
						x = merged,
						file = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-data.tsv'),
						sep = '\t')
					
					data.table::fwrite(
						x = ld,
						file = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-ld.tsv'),
						sep ='\t')
				}
			}
		}
	})

	output$batch_download = downloadHandler(
		filename=function(){return('batch_results.zip')},
		content=function(file){
			if (valid_batch_input()){
				coordinate_list = str_split(trimws(input$batch_input),'\n')[[1]]
			} else {
				coordinate_list = unlist(fread(input$batch_region$datapath,header=FALSE,sep='\t'))
			}
			owd=setwd(tmp_dir)
			on.exit({setwd(owd)})
			tar(file,coordinate_list,compression='gzip')
		}
	)
})
