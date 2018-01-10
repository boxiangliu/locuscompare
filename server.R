# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01

library(shiny)
library(DT)
library(stringr)
library(dplyr)

#----------------------- User Input ---------------------------#


in_fn1='data/eqtl.txt' # input file 1
in_fn2='data/gwas.txt' # input file 2


#--------------------- End of User Input ---------------------#
source('locuscompare.R')

# Variables: 
title1='Study 1' # Study 1 name
title2='Study 2' # Study 2 name
panel='data/integrated_call_samples_v3.20130502.ALL.panel'
tmp_dir='tmp/' # temporary directory
if (!dir.exists(tmp_dir)){dir.create(tmp_dir,recursive=TRUE)}

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

.anno=fread('data/gencode.v19.genes.v6p.hg19.bed',select=c(1:3,5,6),col.names=c('chr','start','end','gene_id','gene_name'))
.anno[,chr:=str_replace(chr,'chr','')]
gene_name2coordinate=function(gene_name_query){
    tmp=.anno[gene_name==gene_name_query,list(chr,start)]
    coordinate=paste0(tmp$chr,':',tmp$start-1e6,'-',tmp$start+1e6)
    return(coordinate)
}

query_database=function(conn,table_name,coordinate='',gene_name=''){
    # validate(need(coordinate!=''|gene_name!=''),'Please enter a coordinate or a gene name')
    # validate(need(coordinate==''|gene_name==''),'Please enter either one but not both')
    if (str_detect(table_name,'^eQTL')) {
        validate(need(gene_name!='','Please enter a gene name for the eQTL dataset'))
    }
    if (gene_name!=''){
        # validate(need(coordinate==''),'Please enter either one but not both')
        coordinate=gene_name2coordinate(gene_name)
    }
    parsed_coordinate=parse_coordinate(coordinate)
    table=tbl(conn,table_name)
    res=table%>%dplyr::filter(chr==parsed_coordinate$chr,pos>=parsed_coordinate$start,pos<=parsed_coordinate$end)%>%collect()
    return(res)
}

shinyServer(function(input, output, session) {    
    print(isolate(input$study1))
    output$debugger=renderText({(input$gene_name=='')})
    
    reactive({
        input$visualize
        # d1=fread(in_fn1)
        # d2=fread(in_fn2)
        d1=isolate(
            query_database(
                conn = locuscompare_db,
                table_name = input$study1,
                coordinate = input$coordinate,
                gene_name = input$gene_name))
        print(isolate(d1))
        d1=isolate(
            query_database(
                conn = locuscompare_db,
                table_name = input$study2,
                coordinate = input$coordinate,
                gene_name = input$gene_name))
        merged=merge(d1,d2,by=c('rsid','chr','pos'),suffixes=c('1','2'),all=FALSE)
        merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
        # retrieve_vcf(merged,tmp_dir)
        chr=unique(merged$chr)
        snp_init=merged[which.min(pval1+pval2),rsid]
        ld_init=NULL
    })

    
    


    # values=reactiveValues(snp_plot=snp_init,
    #                       snp_table=snp_init,
    #                       color=assign_color(merged$rsid,snp_init,ld_init),
    #                       shape=assign_shape(merged,snp_init),
    #                       size=assign_size(merged,snp_init),
    #                       ld=ld_init)

    values=reactiveValues(snp_plot=NULL,
                          snp_table=NULL,
                          color=NULL,
                          shape=NULL,
                          size=NULL,
                          ld=NULL)
    
    observeEvent(input$plot_dblclick,{
        values$snp_plot=select_snp(input$plot_dblclick,merged)
        values$color=assign_color(merged$rsid,values$snp_plot,values$ld)
        values$shape=assign_shape(merged,values$snp_plot)
        values$size=assign_size(merged,values$snp_plot)
    })
    
    observeEvent(input$population,{
        # values$ld=calc_LD(merged$rsid,chr,input$population,tmp_dir,paste0(tmp_dir,'/1000genomes.vcf.gz'),panel)
        values$ld=NULL
        # values$snp_plot=snp_init
        # values$color=assign_color(merged$rsid,snp_init,values$ld)
        # values$shape=assign_shape(merged,snp_init)
        # values$size=assign_size(merged,snp_init)
        values$snp_plot=NULL
        values$color=NULL
        values$shape=NULL
        values$size=NULL
    })
    
    observeEvent(input$plot_click,{
        values$snp_table=select_snp(input$plot_click,merged)
    })
    

    output$locuscompare = renderPlot({
        merged[,label:=ifelse(rsid==values$snp_plot,rsid,'')]
        p1=make_locuscatter(merged,title1,title2,values$ld,values$color,values$shape,values$size,legend=FALSE)
        p1
    },height=function(){
        session$clientData$output_locuscompare_width
    })
    
    output$locuszoom1 = renderPlot({
      p2=make_locuszoom(merged[,list(rsid,logp1,label)],title1,values$ld,values$color,values$shape,values$size,y_string='logp1')
      p2
    })
    
    output$locuszoom2 = renderPlot({
      p3=make_locuszoom(merged[,list(rsid,logp2,label)],title2,values$ld,values$color,values$shape,values$size,y_string='logp2')
      p3
    })
    
    
    output$info = DT::renderDataTable({
      snp_table=values$snp_table
      ld_snps=values$ld[SNP_A==snp_table,][R2>=input$r2_threshold,list(rsid=SNP_B,r2=R2)]
      ld_snps=rbind(data.table(rsid=snp_table,r2=1),ld_snps)
      infoTable=merged[rsid%in%ld_snps$rsid,list(rsid,chr,pos,pval1,pval2)]
      infoTable=merge(infoTable,ld_snps,by='rsid')
      DT::datatable(infoTable)
    })

  
})
