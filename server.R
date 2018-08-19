# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01

options(shiny.maxRequestSize=100*1024^2) 

# Functions:
select_snp=function(click,merged){
	merged_subset=nearPoints(merged,click,maxpoints=1)
	snp=merged_subset %>% dplyr::select(rsid) %>% unlist()
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
			end = chr_start_end$start + gene_window
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
				"select display_trait 
				from %s_trait;",study)
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
	conn <- do.call(DBI::dbConnect, args)
	on.exit(DBI::dbDisconnect(conn))
	if (str_detect(study,'^eQTL')){
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
    
	if (str_detect(study,'^GWAS')){
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
	
	if (valid_study){
		res=dbGetQuery(
			conn = conn,
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
			conn = conn,
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
		res = res %>% dplyr::filter(rsid %in% rsid_list$rsid)
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
		res=fread(datapath,header=TRUE,colClasses=c(trait='character',rsid='character',pval='numeric'))
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

epochTime = function() {
	as.integer(Sys.time())
}

humanTime <- function() {
	format(Sys.time(), "%Y%m%d-%H%M%OS")
}

saveData <- function(data,dir,name) {
	fileName <- sprintf("%s_%s_info.csv",
		name,
		humanTime())

	write.csv(x = data, file = file.path(dir, fileName),
		row.names = FALSE, quote = FALSE)
}

batch_query = function(tmp_dir,coordinate_list,valid_batch_study1,valid_batch_study2,input,token){
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
					plot = p1,
					base_height = input$batch_locuscompare_length,
					base_width = input$batch_locuscompare_length)
				
				cowplot::save_plot(
					filename = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-locuszoom1.pdf'),
					plot = p2,
					base_height = input$batch_locuszoom_height,
					base_width = input$batch_locuszoom_width)
				
				cowplot::save_plot(
					filename = paste0(tmp_dir,'/',coordinate,'/',trait1,'-',trait2,'-locuszoom2.pdf'),
					plot = p3,
					base_height = input$batch_locuszoom_height,
					base_width = input$batch_locuszoom_width)
				
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
	
	tar_fn = paste0(tmp_dir,'/',input$batch_job_name,'-',token,'-',format(Sys.time(), "%Y_%b_%d_%X"),'.tar.gz')

	owd = setwd(tmp_dir)
	suppressWarnings(tar(tar_fn,coordinate_list,compression='gzip'))
	setwd(owd)
	
	return(tar_fn)
}


fwrite_return = function(x,file,sep){
    fwrite(x=x,file=file,sep=sep)
    return(file)
}

ggsave_return = function(filename,plot,width,height){
    ggsave(filename=filename,plot=plot,width=width,height=height)
    return(filename)
}

shinyServer(function(input, output, session) {
	# Session-specific variables:
	tmp_dir = paste0(tempdir(),'/',session$token,'/')
	dir.create(tmp_dir,recursive = TRUE)
	Sys.chmod(tmp_dir, mode="0777")
	hide(id = "loading-content", anim = TRUE, animType = "fade")    
	show("app-content")
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
	    shinyjs::reset('Plots')
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
				end = chr_start_end$start + input$gene_window*1e3
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
		valid_study1_ = valid_study1()
		input_study1_ = input$study1
		input_trait1_ = input$trait1
		input_file1_datapath_ = input$file1$datapath
		coordinate_ = coordinate()
		future({get_study(valid_study1_,input_study1_,input_trait1_,input_file1_datapath_,coordinate_)})
	})

	d2 = eventReactive(input$visualize,{
		shiny::validate(need(valid_study2() | valid_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_study2() & valid_file2()),'Please select or upload study 2, but not both!'))
		valid_study2_ = valid_study2()
		input_study2_ = input$study2
		input_trait2_ = input$trait2
		input_file2_datapath_ = input$file2$datapath
		coordinate_ = coordinate()
		future({get_study(valid_study2_,input_study2_,input_trait2_,input_file2_datapath_,coordinate_)})
	})

	merged=reactive({
	    d1_non_empty = d1() %...>% nrow() %...>% `>`(0)
	    d2_non_empty = d2() %...>% nrow() %...>% `>`(0)
	    shiny::validate(need(d1_non_empty,'No SNP was found in specified region for study 1. Did you input the correct region?'))
	    shiny::validate(need(d1_non_empty,'No SNP was found in specified region for study 2. Did you input the correct region?'))
		merged= promise_all(d1 = d1(), d2= d2()) %...>% {merge(.$d1,.$d2,by='rsid',suffixes=c('1','2'),all=FALSE)}
		merged= merged %...>% get_position()
		check_overlap = merged %...>% nrow() %...>% `>`(0)
		shiny::validate(need(check_overlap,'No overlapping SNPs between two studies'))
		merged = merged %...>% setDT()
		merged = merged %...>% mutate(logp1 = -log10(pval1),logp2 = -log10(pval2))
		return(merged)
	})

	snp=reactiveVal(value='',label='snp')
	
	observeEvent(merged(),{
		# updateSelectizeInput(session, "snp", choices = merged()$rsid, server = TRUE)
		merged() %...>% 
			select(rsid) %...>% 
			unname () %...>%
			unlist() %...>% 
			updateSelectizeInput(session, "snp", choices = ., server = TRUE)
	})
	
	observeEvent(input$visualize,{
	    non_empty_merge = merged() %...>% nrow() %...>% `>`(0)
	    non_empty_merge %...>% (
	        function(non_empty_merge){
	            if (non_empty_merge){
	                merged() %...>% 
	                    dplyr::slice(which.min(pval1*pval2)) %...>% 
	                    dplyr::select(rsid) %...>% 
	                    unlist() %...>% 
	                    snp()
	            } else {
	                snp('')
	            }
	        }
	    )
	})

	observeEvent(input$snp,{
		snp(input$snp)
	})

	chr=reactive({
		chr=unique(coordinate()$chr)
		shiny::validate(need(length(chr)==1,'Studies must only have one chromosome!'))
		return(chr)
	})

	ld=reactive({
		retrieve_LD(chr(),snp(),input$population)
	})

	color=reactive({
		merged() %...>% dplyr::select(rsid) %...>% assign_color(snp(),ld())
	})
	shape=reactive({
		merged() %...>% assign_shape(snp())
	})
	size=reactive({
		merged() %...>% assign_size(snp())
	})

	observeEvent(input$plot_click,{
		selected_snp = merged() %...>% select_snp(input$plot_click,.)
		nonempty_snp = selected_snp %...>% identical(character(0)) %...>% `!`
		nonempty_snp %...>% (
			function(nonempty_snp){
				if (nonempty_snp){
					selected_snp %...>% snp() 
				}
			}
		)
	})

	range=reactiveValues(xmin=NULL,xmax=NULL)

	observeEvent(input$plot_brush,{
	    brush = input$plot_brush
        range$xmin=brush$xmin
        range$xmax=brush$xmax
	})
	
	observeEvent(input$plot_dblclick,{
        range$xmin=NULL
        range$xmax=NULL
	})
	
	plot_data=reactive({
		if (is.null(range$xmin)){
			plot_data = merged()
		} else {
			plot_data = merged() %...>% dplyr::filter(pos<=range$xmax,pos>=range$xmin)
		}
		plot_data = plot_data %...>% mutate(label=ifelse(rsid==snp(),rsid,''))
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
		p = promise_all(plot_data = plot_data(), color = color(),shape = shape(), size = size()) %...>% {
			make_locuscatter(
			merged = .$plot_data,
			title1 = title1(),
			title2 = title2(),
			ld = ld(),
			color = .$color,
			shape = .$shape,
			size = .$size,
			legend=FALSE
			)
		}
		return(p)
	})
	
	output$locuszoom1 = renderPlot({
		p = promise_all(plot_data = plot_data(), color = color(),shape = shape(), size = size()) %...>% {
			make_locuszoom(
			metal = .$plot_data,
			title = title1(),
			ld = ld(),
			color = .$color,
			shape = .$shape,
			size = .$size,
			y_string='logp1')
		}
		return(p)
	})
	
	output$locuszoom2 = renderPlot({
		p = promise_all(plot_data = plot_data(), color = color(),shape = shape(), size = size()) %...>% {
			make_locuszoom(
			metal = .$plot_data,
			title = title2(),
			ld = ld(),
			color = .$color,
			shape = .$shape,
			size = .$size,
			y_string='logp2')
		}
		return(p)
	})
	
	output$blank_plot = renderPlot({
		p = plot_data() %...>% {ggplot() + geom_blank()}
		return(p)
	})
	
	output$snp_info = renderText({
		res = dbGetQuery(
			conn = locuscompare_pool,
			statement = sprintf('select * from tkg_p3v5a where rsid = "%s";',snp()))
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
	
	# TODO: Use this code when DT is compatible with async.
	# output$ld_snps = DT::renderDataTable({
	# 	ld_snps=ld() %>%
	# 		dplyr::filter(SNP_A==snp(),R2>=input$r2_threshold) %>%
	# 		dplyr::select(rsid=SNP_B,r2=R2)
	# 	ld_snps=rbind(data.frame(rsid=snp(),r2=1),ld_snps)
	# 	ld_snps_2 = merged() %...>%
	# 		dplyr::select(rsid,chr,pos,pval1,pval2) %...>%
	# 		merge(ld_snps,by='rsid')
	#   return(ld_snps_2)
	# })

	output$ld_snps = renderTable({
		ld_snps=ld() %>%
			dplyr::filter(SNP_A==snp(),R2>=input$r2_threshold) %>%
			dplyr::select(rsid=SNP_B,r2=R2)
		ld_snps=rbind(data.frame(rsid=snp(),r2=1),ld_snps)
		ld_snps_2 = merged() %...>%
		    dplyr::mutate(pval1_disp = format.pval(pval1), pval2_disp = format.pval(pval2)) %...>%
			dplyr::select(rsid,chr,pos,pval1_disp,pval2_disp) %...>%
			merge(ld_snps,by='rsid') %...>%
		    dplyr::rename(rsID = rsid, Chromosome = chr, Position = pos, `P-value 1` = pval1_disp, `P-value 2` = pval2_disp)
		return(ld_snps_2)
	},width = '100%', striped = TRUE, hover = TRUE, bordered = TRUE)
	
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
			data_fn = merged() %...>% fwrite_return(paste0(tmp_dir,'/data.tsv'),sep='\t')
			fwrite(ld(),paste0(tmp_dir,'/ld.tsv'),sep='\t')

			locuscompare = promise_all(plot_data = plot_data(), color = color(),shape = shape(), size = size()) %...>% {
			    make_locuscatter(
			        merged = .$plot_data,
			        title1 = title1(),
			        title2 = title2(),
			        ld = ld(),
			        color = .$color,
			        shape = .$shape,
			        size = .$size,
			        legend=FALSE)}
			
			locuszoom1 = promise_all(plot_data = plot_data(), color = color(),shape = shape(), size = size()) %...>% {
			    make_locuszoom(
			        metal = .$plot_data,
			        title = title1(),
			        ld = ld(),
			        color = .$color,
			        shape = .$shape,
			        size = .$size,
			        y_string='logp1')}
			    
			locuszoom2 = promise_all(plot_data = plot_data(), color = color(),shape = shape(), size = size()) %...>% {
			    make_locuszoom(
			        metal = .$plot_data,
			        title = title2(),
			        ld = ld(),
			        color = .$color,
			        shape = .$shape,
			        size = .$size,
			        y_string='logp2')}

			
			length_ = input$locuscompare_length
			width_ = input$locuszoom_width
			height_ = input$locuszoom_height
			locuscompare_fn = locuscompare %...>% ggsave_return(paste0(tmp_dir,'/locuscompare.jpg'),.,width=length_,height=length_)
			locuszoom1_fn = locuszoom1 %...>% ggsave_return(paste0(tmp_dir,'/locuszoom1.jpg'),.,width=width_,height=height_)
			locuszoom2_fn = locuszoom2 %...>% ggsave_return(paste0(tmp_dir,'/locuszoom2.jpg'),.,width=width_,height=height_)
			browser()
			promise_all(data_fn = data_fn, locuscompare_fn = locuscompare_fn, locuszoom1_fn = locuszoom1_fn, locuszoom2_fn = locuszoom2_fn) %...>% {
			    c(.$data_fn,paste0(tmp_dir,'/ld.tsv'),.$locuscompare_fn,.$locuszoom1_fn,.$locuszoom2_fn) %>% utils::zip(file,.,flags = '-j')}
		}
		
	)

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
        
		shiny::validate(need(isTruthy(input$batch_job_name),'Please provide job name!'))
		shiny::validate(need(isTruthy(input$batch_job_email),'Please provide an email!'))

		return('All inputs are given. Ready to submit!')
	})

	observeEvent(input$submit,{
		shiny::validate(need(valid_batch_study1() | valid_batch_file1(),'Please provide study 1!'))
		shiny::validate(need(!(valid_batch_study1() & valid_batch_file1()),'Please select or upload study 1, but not both!'))
		
		shiny::validate(need(valid_batch_study2() | valid_batch_file2(),'Please provide study 2!'))
		shiny::validate(need(!(valid_batch_study2() & valid_batch_file2()),'Please select or upload study 2, but not both!'))

		shiny::validate(need(any(valid_batch_input(),valid_batch_region()),'Please provide a list of regions!'))
		shiny::validate(need({valid_batch_input()+valid_batch_region()==1},'Please either input or upload a list of regions, but not both!'))
		
		if (valid_batch_input()){
			coordinate_list = str_split(trimws(input$batch_input),'\n')[[1]]
		} else {
			coordinate_list = unlist(fread(input$batch_region$datapath,header=FALSE,sep='\t'))
		}

		coordinate_list = coordinate_list[1:min(25,length(coordinate_list))] # Subset to first 25 because of gmail size limit.
		
		valid_batch_study1_ = function() {valid_batch_study1()}
		valid_batch_study2_ = function() {valid_batch_study2()}
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
	
	observeEvent(input$submit,{
		shinyjs::hide('batch_query')
		shinyjs::show('batch_query_success')
	})
	
	
	observeEvent(input$submit_another_query,{
		shinyjs::reset('batch_query')
		shinyjs::show('batch_query')
		shinyjs::hide('batch_query_success')
	})
	
	#---------------# 
	# Download page #
	#---------------#
	sheet_key = googlesheets::gs_key(x='1gq46xlOk674Li50cpv9riZYG7bsfSeZB5qSefa82bR8',lookup=FALSE)
	list_of_studies = googlesheets::gs_read(sheet_key)
	output$study_info = renderDataTable({
		DT::datatable(list_of_studies)
	})

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
