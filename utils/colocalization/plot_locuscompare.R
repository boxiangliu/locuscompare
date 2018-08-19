# detach('package:locuscomparer',unload=TRUE)
# devtools::install_github("boxiangliu/locuscomparer")
library(data.table)
library(locuscomparer)
source('global.R')
source('utils/utils/genome_annotation.R')
library(ggrepel)
library(cowplot)

fig_dir = 'figure/colocalization/plot_locuscompare/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
out_dir = 'processed_data/colocalization/plot_locuscompare/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

gencode = read_gencode()

gwas_tbl = 'GWAS_Coronary_Artery_Disease_Nelson_2017'
gwas_trait = 'Coronary-Artery-Disease'
eqtl_tbl = 'eQTL_Artery_Coronary_GTEx_v6p'
gwas_title = 'CAD GWAS'
eqtl_title = 'Coronary Artery eQTL'

both = c('MRAS','PHACTR1')
ecaviar = c('SEMA5A','FES','ARVCF')
smr = c('NBEAL1','SYPL2','ADAMTS7')


get_gwas = function(gwas_tbl,gwas_trait,gene_coordinate){
	gwas = dbGetQuery(
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
			gwas_tbl,
			gwas_trait,
			gene_coordinate$chr,
			gene_coordinate$start,
			gene_coordinate$end
			)
		)
	setDT(gwas)
	return(gwas)
}

get_eqtl = function(eqtl_tbl,gene_id,gene_coordinate){
	eqtl = dbGetQuery(
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
			eqtl_tbl,
			gene_id,
			gene_coordinate$chr,
			gene_coordinate$start,
			gene_coordinate$end
			)
		)
	setDT(eqtl)
	return(eqtl)
}

make_plot = function(gene_set,fig_dir,return=FALSE,top_snp=NULL){
	p_list = list()
	for(i in seq_along(gene_set)){
		gene = gencode[gene_name == gene_set[i],]
		gene_coordinate = list(
			chr = str_replace(gene$chr,'chr',''),
			start = gene$start - 1e6,
			end = gene$start + 1e6)
		gene_id = gene$gene_id 
		gene_name = gene$gene_name
		message(gene_name)

		gwas = get_gwas(gwas_tbl,gwas_trait,gene_coordinate)
		eqtl = get_eqtl(eqtl_tbl,gene_id,gene_coordinate)
		merged = merge(gwas,eqtl,by='rsid')
		merged[,prod:=pval.x*pval.y]

		if (is.na(top_snp[i])){
			lead_SNP = merged[which.min(merged$prod),rsid]
		} else {
			lead_SNP = top_snp[i]
		}
		
		vcf_fn = sprintf('%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir,gene_coordinate$chr)
		chr_start_end = sprintf('%s:%s-%s',gene_coordinate$chr,gene_coordinate$start,gene_coordinate$end)
		vcf_out_fn = sprintf('%s/%s.vcf.gz',out_dir,chr_start_end)
		command = sprintf('bcftools view %s %s -Oz -o %s',vcf_fn,chr_start_end,vcf_out_fn)
		system(command)

		p1 = locuscomparer::main(
			in_fn1 = gwas,
			in_fn2 = eqtl,
			title1 = gwas_title,
			title2 = eqtl_title,
			vcf_fn = vcf_out_fn,
			combine = TRUE,
			snp = lead_SNP,
			legend_position = 'bottomright',
			lz_ylab_linebreak = TRUE
			)

		if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
		fig_fn = sprintf('%s/%s_%s.pdf',fig_dir,gene_name,lead_SNP)
		save_plot(fig_fn,p1,base_height=4,base_width=8)
		p_list[[gene_name]] = ggdraw(p1) + draw_label(gene_name,x = 0.2, y = 0.95)
	}
	if (return){
		return(p_list)
	}
}

p_list_1 = make_plot(both,paste0(fig_dir,'/both'),return=TRUE)
p_grid_1 = plot_grid(plotlist = p_list_1, ncol=1, labels = letters[1:length(p_list_1)])
fig_fn = sprintf('%s/both/both_grid.pdf',fig_dir)
save_plot(fig_fn,p_grid_1,base_width=8,base_height=4*length(p_list_1))


p_list_2 = make_plot(ecaviar,paste0(fig_dir,'/ecaviar'),return=TRUE)
p_grid_2 = plot_grid(plotlist = p_list_2, ncol=1, labels = letters[1:length(p_list_2)])
fig_fn = sprintf('%s/ecaviar/ecaviar_grid.pdf',fig_dir)
save_plot(fig_fn,p_grid_2,base_width=8,base_height=4*length(p_list_2))

p_list_3 = make_plot(smr,paste0(fig_dir,'/smr'),return=TRUE,top_snp=c(NA,'rs2272272',NA))
p_grid_3 = plot_grid(plotlist = p_list_3, ncol=1, labels = letters[1:length(p_list_3)])
fig_fn = sprintf('%s/smr/smr_grid.pdf',fig_dir)
save_plot(fig_fn,p_grid_3,base_width=8,base_height=4*length(p_list_3))
