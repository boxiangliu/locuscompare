# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01

library(shiny)
library(DT)
library(stringr)
library(dplyr)
source('locuscompare.R')
library(RMySQL)
home_dir='/srv/persistent/bliu2/locuscompare/'


# Variables: 
title1='Study 1' # Study 1 name
title2='Study 2' # Study 2 name
panel=paste0(home_dir,'data/integrated_call_samples_v3.20130502.ALL.panel')
tmp_dir=paste0(home_dir,'tmp/') # temporary directory
if (!dir.exists(tmp_dir)){dir.create(tmp_dir,recursive=TRUE)}

locuscompare_db <- dbPool(
    RMySQL::MySQL(), 
    dbname = "locuscompare",
    host = "rds-mysql-locuscompare.cbhpzvkzr3rc.us-west-1.rds.amazonaws.com",
    username = "admin",
    password = "12345678"
)

# Functions:
select_snp=function(click,merged){
    merged_subset=nearPoints(merged,click)
    if (nrow(merged_subset)>1) {
        merged_subset=merged_subset[1,]
    }
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

query_database=function(conn,table_name,locus){
    if (str_detect(table_name,'^eQTL')) {
        shiny::validate(need(is_coordinate(locus)==FALSE,'Since an eQTL dataset is selected, you must enter a gene name'))
    }
    if (!is_coordinate(locus)){
        coordinate=gene_name2coordinate(locus)
    }
    parsed_coordinate=parse_coordinate(coordinate)
    table=tbl(conn,table_name)
    res=table%>%dplyr::filter(chr==parsed_coordinate$chr,pos>=parsed_coordinate$start,pos<=parsed_coordinate$end)%>%collect()
    return(res)
}

shinyServer(function(input, output, session) {
    
    output$exampleInput=renderTable({
        data.frame(rsid=c('rs56171612','rs1866475','rs1866476'),
                   chr='15',
                   pos=as.integer(c(90427941,90428072,90428543)),
                   pval=c(0.346735,0.644105,0.797991))
    })
    
    merged=eventReactive(input$visualize,{
        shiny::validate(need(input$study1!='' | !is.null(input$upload_study1),'Please upload or select study 1'))
        shiny::validate(need(input$study2!='' | !is.null(input$upload_study2),'Please upload or select study 2'))
        shiny::req(input$locus)
        
        if (input$study1!=''){
            d1=query_database(
                conn = locuscompare_db,
                table_name = input$study1,
                locus = input$locus)
        } else {
            d1=fread(input$upload_study1$datapath)
        }
    
        if (input$study2!=''){
            d2=query_database(
                conn = locuscompare_db,
                table_name = input$study2,
                locus = input$locus)
        } else {
            d2=fread(input$upload_study2$datapath)
        }

        merged=merge(d1,d2,by=c('rsid','chr','pos'),suffixes=c('1','2'),all=FALSE)
        setDT(merged)
        merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
        return(merged)
    })
    
    snp=reactiveVal(value=NULL,label='snp')
    
    observeEvent(input$visualize,{
        retrieve_vcf(merged(),tmp_dir)
        snp(merged()[which.min(pval1+pval2),rsid])
    })

    chr=reactive({unique(merged()$chr)})
    ld=reactive({calc_LD(merged()$rsid,chr(),input$population,tmp_dir,paste0(tmp_dir,'/1000genomes.vcf.gz'),panel)})
    
    observeEvent(input$plot_click,{
        snp(select_snp(input$plot_click,merged()))
    })

    color=reactive({assign_color(merged()$rsid,snp(),ld())})
    shape=reactive({assign_shape(merged(),snp())})
    size=reactive({assign_size(merged(),snp())})
    
    plot_data=reactive({
        plot_data=merged()
        plot_data[,label:=ifelse(rsid==snp(),rsid,'')]
        return(plot_data)
    })

    locuscompare=reactive({
        make_locuscatter(merged = plot_data(),
                         title1 = input$study1,
                         title2 = input$study2,
                         ld = ld(),
                         color = color(),
                         shape = shape(),
                         size = size(),
                         legend=FALSE)
    })
    
    locuszoom1=reactive({
        make_locuszoom(
            metal = plot_data()[,list(rsid,logp1,label)],
            title = input$study1,
            ld = ld(),
            color = color(),
            shape = shape(),
            size = size(),
            y_string='logp1')
    })

    locuszoom2=reactive({
        make_locuszoom(
            metal = plot_data()[,list(rsid,logp2,label)],
            title = input$study2,
            ld = ld(),
            color = color(),
            shape = shape(),
            size = size(),
            y_string='logp2')
    })
    
    output$locuscompare = renderPlot({
        locuscompare()
    },height=function(){
        session$clientData$output_locuscompare_width
    })
    
    output$locuszoom1 = renderPlot({
        locuszoom1()
    })
    
    output$locuszoom2 = renderPlot({
        locuszoom2()
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
            
            zip(file,c('data.tsv','ld.txt','locuscompare.pdf','locuszoom_study1.pdf','locuszoom_study2.pdf'))
        }
        
    )
})
