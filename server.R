# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01


# Functions:
select_snp=function(click,merged){
	merged_subset=nearPoints(merged,click,maxpoints=1)
	snp=merged_subset[,rsid]
	return(snp)
}

parse_coordinate=function(coordinate){
	split_coordinate=str_split_fixed(coordinate,':',2)
	chr=split_coordinate[,1]
	chr=str_replace(chr,'chr','')
	pos=split_coordinate[,2]
	split_pos=str_split_fixed(pos,'-',2)
	start=as.integer(split_pos[,1])
	end=as.integer(split_pos[,2])
	return(list(chr=chr,start=start,end=end))
}

.anno=fread(paste0(home_dir,'data/gencode/gencode.v19.genes.v6p.hg19.bed'),select=c(1:3,5,6),col.names=c('chr','start','end','gene_id','gene_name'))
.anno[,chr:=str_replace(chr,'chr','')]
gene_name2coordinate=function(gene_name_query){
	tmp=.anno[gene_name==gene_name_query,list(chr,start)]
	coordinate=paste0(tmp$chr,':',tmp$start-1e5,'-',tmp$start+1e5)
	return(coordinate)
}

is_coordinate=function(locus){
	if (str_detect(locus,'^chr|^[0-9]{1}')){
		return(TRUE)
	} else {
		return(FALSE)
	}
}

get_trait=function(study){
	if (study==''){
		trait=''
	} else if (str_detect(study,'^GWAS_')){
		x = paste0(study,'_traits')
		table=tbl(locuscompare_db,x)
		res=table%>%dplyr::select(trait)%>%collect()%>%unlist()%>%unname()
		trait=sort(res)
	} else if (str_detect(study,'^eQTL_')){
		x = paste0(study,'_traits')
		table=tbl(locuscompare_db,x)
		res=table%>%dplyr::select(gene_name)%>%collect()%>%unlist()%>%unname()
		trait=sort(res)
	} else {
		NULL
	}
	return(trait)
}

query_database=function(conn,table_name,locus,trait){
	if (!is_coordinate(locus)){
		locus=gene_name2coordinate(locus)
	}
	locus=str_replace_all(locus,',','')
	parsed_coordinate=parse_coordinate(locus)
	table=tbl(conn,table_name)
	
	if (trait=='batch'){ # for batch mode
		res=table%>%
			dplyr::filter(chr==parsed_coordinate$chr,
						  pos>=parsed_coordinate$start,
						  pos<=parsed_coordinate$end)%>%collect()
		setDT(res)
		if ('gene_name' %in% colnames(res)){
			setnames(res,'gene_name','trait')
		}
		stopifnot('trait' %in% colnames(res))
	} else {
		if (str_detect(table_name,'^eQTL_')) {
			print(table_name)
			res=table%>%
				dplyr::filter(gene_id==trait,
							  chr==parsed_coordinate$chr,
							  pos>=parsed_coordinate$start,
							  pos<=parsed_coordinate$end)%>%collect()
			print('finished reading table')
		} else if (str_detect(table_name,'^GWAS_')){
			res=table%>%
				dplyr::filter(chr==parsed_coordinate$chr,
							  pos>=parsed_coordinate$start,
							  pos<=parsed_coordinate$end)%>%collect()
		}
	}

	return(res)
}

get_study = function(valid_study,study,trait,datapath,coordinate){
	if (valid_study){
		if (str_detect(study,'^eQTL')){
			column='gene_name'
		} else {
			column='trait'
		}
		res=dbGetQuery(
			conn = locuscompare_db,
			statement = sprintf(
				'select rsid,chr,pos,pval from %s where %s = "%s" and chr = "%s" and pos >= %s and pos <= %s;',
				study,
				column,
				trait,
				coordinate$chr,
				coordinate$start,
				coordinate$end
				)
			)
	} else {
		res=fread(datapath,header=TRUE,colClasses=c(rsid='character',chr='character',pos='integer',pval='numeric'))
		shiny::validate(need(all(c('rsid','chr','pos','pval')%in%colnames(res)),'Input file must have columns rsid, chr, pos, pval!'))
	}
	return(res)
}

shinyServer(function(input, output, session) {
	# Update trait field according to selected study:
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
				conn = reference_db,
				statement = sprintf('select chr,pos from tkg_p3v5a where rsid = "%s";',input$reference_snp)
				)
			res=list(
				chr = chr_pos$chr,
				start = chr_pos$pos - input$snp_window*1e3,
				end = chr_pos$pos + input$snp_window*1e3
				)
		}
		if (valid_gene_region()){
			chr_start_end=dbGetQuery(
				conn = reference_db,
				statement = sprintf('select chr,start,end from gencode_v19_gtex_v6p where gene_name = "%s";',input$reference_gene)
				)
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
		merged=merge(d1(),d2(),by=c('rsid','chr','pos'),suffixes=c('1','2'),all=FALSE)
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
		chr=unique(merged()$chr)
		validate(need(length(chr)==1,'Studies must only have one chromosome!'))
	})

	vcf_fn=reactive({
		retrieve_vcf(merged(),tmp_dir)
	})

	ld=reactive({
		calc_LD(merged()$rsid,chr(),input$population,tmp_dir,vcf_fn())
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
	
	output$locuscompare = renderPlot({
		make_locuscatter(
			merged = plot_data(),
			title1 = input$trait1,
			title2 = input$trait2,
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
			title = input$trait1,
			ld = ld(),
			color = color(),
			shape = shape(),
			size = size(),
			y_string='logp1')
	})
	
	output$locuszoom2 = renderPlot({
		make_locuszoom(
			metal = plot_data()[,list(rsid,chr,pos,logp2,label)],
			title = input$trait2,
			ld = ld(),
			color = color(),
			shape = shape(),
			size = size(),
			y_string='logp2')
	})
	
	output$snp_info = renderText({
		res=dbGetQuery(
			conn = reference_db,
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
		filename = function(){return('eQTL_Lung_GTEx_2017_IL18R1.tsv')},
		content = function(file){file.copy('data/example/eQTL_Lung_GTEx_2017_IL18R1.tsv',file)},
		contentType = 'text/txt')

	output$file2_example = downloadHandler(
		filename = function(){return('GWAS_Asthma_Moffatt_2010_chr2_101973286_103971063.tsv')},
		content = function(file){file.copy('data/example/GWAS_Asthma_Moffatt_2010_chr2_101973286_103971063.tsv',file)},
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
				title1 = input$trait1,
				title2 = input$trait2,
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				legend=FALSE
				)

			locuszoom1 = make_locuszoom(
				metal = plot_data()[,list(rsid,chr,pos,logp1,label)],
				title = input$trait1,
				ld = ld(),
				color = color(),
				shape = shape(),
				size = size(),
				y_string='logp1'
				)

			locuszoom2 = make_locuszoom(
				metal = plot_data()[,list(rsid,chr,pos,logp2,label)],
				title = input$trait2,
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

	output$batch_file1_example = downloadHandler(
		filename = function(){return('GWAS_Asthma_Moffatt_2010_chr2_chr9_chr22.tsv')},
		content = function(file){file.copy('data/example/GWAS_Asthma_Moffatt_2010_chr2_chr9_chr22.tsv',file)},
		contentType = 'text/txt'
		)

	output$batch_file2_example = downloadHandler(
		filename = function(){return('eQTL_Lung_GTEx_2017_IL18R1_IL33_IL2RB.tsv')},
		content = function(file){file.copy('data/example/eQTL_Lung_GTEx_2017_IL18R1_IL33_IL2RB.tsv',file)},
		contentType = 'text/txt'
		)

	output$batch_region_example = downloadHandler(
		filename = function(){return('batch_region.txt')},
		content = function(file){file.copy('data/example/batch_region.txt',file)},
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
	})

	observeEvent(input$submit,{
		split_coordinate_input=str_split(trimws(input$batch_coordinate_input),'\n')[[1]]
		
		# Create a Progress object
		progress = shiny::Progress$new()
		on.exit(progress$close())
		progress$set(message = "Making plot", value = 0)
		n = length(split_coordinate_input)
		warning(n)
		for (coordinate in split_coordinate_input){
			warning(tmp_dir)
			warning(coordinate)
			dir.create(paste0(tmp_dir,'/',coordinate),recursive=TRUE)
			parsed_coordinate=parse_coordinate(coordinate)
			
			# Increment the progress bar, and update the detail text.
			progress$inc(1/n, detail = paste("coordinate: ", coordinate))


			d1=query_database(
				conn = locuscompare_db,
				table_name = input$study1,
				locus = coordinate,
				trait = 'batch')
			d2=query_database(
				conn = locuscompare_db,
				table_name = input$study2,
				locus = coordinate,
				trait = 'batch')
			
			trait1_list=unique(d1$trait)
			trait2_list=unique(d2$trait)
			
			for (trait1 in trait1_list){
				for (trait2 in trait2_list){
					print(paste0('trait1: ',trait1,'; trait2: ',trait2))
					merged=merge(d1[trait==trait1],d2[trait==trait2,],by=c('rsid','chr','pos'),suffixes=c('1','2'),all=FALSE)
					if (nrow(merged)==0) next # skip if no overlap
					setDT(merged)
					merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
					
					vcf_fn=retrieve_vcf(merged,tmp_dir)
					snp=merged[which.min(pval1*pval2),rsid]
					chr=parsed_coordinate$chr
					ld=calc_LD(merged$rsid,chr,input$population,tmp_dir,vcf_fn)
					color=assign_color(merged$rsid,snp,ld)
					shape=assign_shape(merged,snp)
					size=assign_size(merged,snp)
					
					plot_data=merged
					plot_data[,label:=ifelse(rsid==snp,rsid,'')]
					
					p1=make_locuscatter(
						merged = plot_data,
						title1 = input$trait1,
						title2 = input$trait2,
						ld = ld,
						color = color,
						shape = shape,
						size = size,
						legend=FALSE)
					
					p2=make_locuszoom(
						metal = plot_data[,list(rsid,chr,pos,logp1,label)],
						title = input$trait1,
						ld = ld,
						color = color,
						shape = shape,
						size = size,
						y_string='logp1')
					
					
					p3=make_locuszoom(
						metal = plot_data[,list(rsid,chr,pos,logp2,label)],
						title = input$trait2,
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
			owd=setwd(tmp_dir)
			on.exit({setwd(owd)})
			split_coordinate_input=str_split(input$batch_coordinate_input,'\n')[[1]]
			zip::zip(file,split_coordinate_input)
		}
	)
	

})
