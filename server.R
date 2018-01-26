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

.anno=fread(paste0(home_dir,'data/gencode.v19.genes.v6p.hg19.bed'),select=c(1:3,5,6),col.names=c('chr','start','end','gene_id','gene_name'))
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

get_traits=function(conn,table_name){

  return(res)
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
    res=table%>%dplyr::select(gene_id)%>%collect()%>%unlist()%>%unname()
    trait=sort(res)
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

shinyServer(function(input, output, session) {
    
    output$exampleInput=renderTable({
        data.frame(rsid=c('rs56171612','rs1866475','rs1866476'),
                   chr='15',
                   pos=as.integer(c(90427941,90428072,90428543)),
                   pval=c(0.346735,0.644105,0.797991))
    })
    
    observe({
      study1_trait=get_trait(input$study1)
      updateSelectizeInput(session, "study1_trait", choices = study1_trait, server = TRUE)
    })
    
    observe({
      study2_trait=get_trait(input$study2)
      updateSelectizeInput(session, "study2_trait", choices = study2_trait, server = TRUE)
    })
    
    merged=eventReactive(input$visualize,{
        shiny::validate(need(input$study1!='' | !is.null(input$upload_study1),'Please upload or select study 1'))
        shiny::validate(need(input$study2!='' | !is.null(input$upload_study2),'Please upload or select study 2'))
        shiny::req(input$locus)
        if (input$study1!=''){
            d1=query_database(
                conn = locuscompare_db,
                table_name = input$study1,
                locus = input$locus,
                trait = input$study1_trait)
        } else {
            d1=fread(input$upload_study1$datapath)
        }
    
        if (input$study2!=''){
            d2=query_database(
                conn = locuscompare_db,
                table_name = input$study2,
                locus = input$locus,
                trait = input$study2_trait)
        } else {
            d2=fread(input$upload_study2$datapath)
        }

        merged=merge(d1,d2,by=c('rsid','chr','pos'),suffixes=c('1','2'),all=FALSE)
        shiny::validate(need(nrow(merged)>0,'No overlapping SNPs between two studies'))
        setDT(merged)
        merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
        return(merged)
    })
    
    
    
    vcf_fn=reactive({retrieve_vcf(merged(),tmp_dir)})
    
    snp=reactiveVal(value=NULL,label='snp')
    
    observeEvent(input$visualize,{
        snp(merged()[which.min(pval1*pval2),rsid])
    })

    chr=reactive({unique(merged()$chr)})
    ld=reactive({calc_LD(merged()$rsid,chr(),input$population,tmp_dir,vcf_fn())})
    
    observeEvent(input$plot_click,{
        selected_snp=select_snp(input$plot_click,merged())
        if (!identical(selected_snp,character(0))){
            snp(selected_snp)    
        }
    })

    color=reactive({assign_color(merged()$rsid,snp(),ld())})
    shape=reactive({assign_shape(merged(),snp())})
    size=reactive({assign_size(merged(),snp())})
    
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
        make_locuscatter(merged = plot_data(),
                         title1 = input$study1,
                         title2 = input$study2,
                         ld = ld(),
                         color = color(),
                         shape = shape(),
                         size = size(),
                         legend=FALSE)
    },height=function(){
        session$clientData$output_locuscompare_width
    })
    
    output$locuszoom1 = renderPlot({
        make_locuszoom(
            metal = plot_data()[,list(rsid,chr,pos,logp1,label)],
            title = input$study1,
            ld = ld(),
            color = color(),
            shape = shape(),
            size = size(),
            y_string='logp1')
    })
    
    output$locuszoom2 = renderPlot({
        make_locuszoom(
            metal = plot_data()[,list(rsid,chr,pos,logp2,label)],
            title = input$study2,
            ld = ld(),
            color = color(),
            shape = shape(),
            size = size(),
            y_string='logp2')
    })
    
    
    output$snp_info = DT::renderDataTable({
      shiny::validate(
          need(input$r2_threshold>0.2,'r2 threshold must be greater than 0.2')
      )
      ld_snps=ld()[SNP_A==snp(),][R2>=input$r2_threshold,list(rsid=SNP_B,r2=R2)]
      ld_snps=rbind(data.table(rsid=snp(),r2=1),ld_snps)
      tmp=merged()[rsid%in%ld_snps$rsid,list(rsid,chr,pos,pval1,pval2)]
      snp_info=merge(tmp,ld_snps,by='rsid')
      DT::datatable(snp_info)
    })
    
    output$download = downloadHandler(
        filename=function(){return('results.zip')},
        content=function(file){
            owd=setwd(tempdir<-tempdir())
            on.exit({setwd(owd)})
            
            fwrite(merged(),'data.tsv',sep='\t')
            fwrite(ld(),'ld.txt',sep='\t')
            ggsave('locuscompare.pdf',locuscompare())
            ggsave('locuszoom_study1.pdf',locuszoom1())
            ggsave('locuszoom_study2.pdf',locuszoom2())
            
            zip(file,c('data.tsv','ld.tsv','locuscompare.pdf','locuszoom_study1.pdf','locuszoom_study2.pdf'))
        }
        
    )


    observeEvent(input$submit_batch_coordinate,{
        print(input$batch_coordinate_input)
        warning(trimws(input$batch_coordinate_input))
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
                        title1 = input$study1,
                        title2 = input$study2,
                        ld = ld,
                        color = color,
                        shape = shape,
                        size = size,
                        legend=FALSE)
                    
                    p2=make_locuszoom(
                        metal = plot_data[,list(rsid,chr,pos,logp1,label)],
                        title = input$study1,
                        ld = ld,
                        color = color,
                        shape = shape,
                        size = size,
                        y_string='logp1')
                    
                    
                    p3=make_locuszoom(
                        metal = plot_data[,list(rsid,chr,pos,logp2,label)],
                        title = input$study2,
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
