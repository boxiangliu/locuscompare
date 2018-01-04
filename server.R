# Server file for LocusCompare
# Boxiang Liu
# 2018-01-01

library(shiny)
library(DT)
source('locuscompare.R')


#----------------------- User Input ---------------------------#
tmp_dir='tmp/' # temporary directory
in_fn1='data/eqtl.txt' # input file 1
in_fn2='data/gwas.txt' # input file 2
tabix='/Users/boshliu/Documents/tools/htslib/bin/tabix' # path to tabix
bgzip='/Users/boshliu/Documents/tools/htslib/bin/bgzip' # path to bgzip 
plink='/Users/boshliu/Documents/tools/plink_mac//plink' # path to plink 
title1='Study 1' # Study 1 name
title2='Study 2' # Study 2 name

#--------------------- End of User Input ---------------------#


# Variables: 
panel='data/integrated_call_samples_v3.20130502.ALL.panel' 
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


shinyServer(function(input, output, session) {    

    d1=fread(in_fn1)
    d2=fread(in_fn2)
    merged=merge(d1,d2,by=c('rsid','chr','pos'),suffixes=c('1','2'),all=FALSE)
    merged[,c('logp1','logp2'):=list(-log10(pval1),-log10(pval2))]
    retrieve_vcf(merged,tmp_dir)

    chr=unique(merged$chr)
    
    snp_init=merged[which.min(pval1+pval2),rsid]
    ld_init=NULL

    values=reactiveValues(snp_plot=snp_init,
                          snp_table=snp_init,
                          color=assign_color(merged$rsid,snp_init,ld_init),
                          shape=assign_shape(merged,snp_init),
                          size=assign_size(merged,snp_init),
                          ld=ld_init)

    observeEvent(input$plot_dblclick,{
        values$snp_plot=select_snp(input$plot_dblclick,merged)
        values$color=assign_color(merged$rsid,values$snp_plot,values$ld)
        values$shape=assign_shape(merged,values$snp_plot)
        values$size=assign_size(merged,values$snp_plot)
    })
    
    observeEvent(input$population,{
        values$ld=calc_LD(merged$rsid,chr,input$population,tmp_dir,paste0(tmp_dir,'/1000genomes.vcf.gz'),panel)
        values$snp_plot=snp_init
        values$color=assign_color(merged$rsid,snp_init,values$ld)
        values$shape=assign_shape(merged,snp_init)
        values$size=assign_size(merged,snp_init)
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
