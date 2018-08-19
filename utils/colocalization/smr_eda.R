library(data.table)
library(stringr)
library(cowplot)
source('utils/utils/genome_annotation.R')

smr_fn = 'processed_data/colocalization/Nelson_gtexv7_coronaryartery_pvalthreshold5e-5.out.smr'
ukbb_fn = '/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'

gencode = read_gencode()
smr = fread(smr_fn)
smr$gene_type = gencode$gene_type[match(smr$Gene,gencode$gene_name)]
smr = smr[gene_type == 'protein_coding']
smr[, padj_SMR := p.adjust(p_SMR,method='fdr')]
setorder(smr,p_SMR)
smr[padj_SMR < 0.05] # nominal p-value cutoff = 3.390201e-05


sypl2 = fread('utils/colocalization/SYPL2.txt')
sypl2[,pos:=as.numeric(str_split_fixed(V2,'_',5)[,2])]
ggplot(sypl2,aes(pos,-log10(V7)))+geom_point()

gwas = fread(sprintf('zcat %s',ukbb_fn))
merged = merge(sypl2,gwas,by.x='pos',by.y='bp_hg19')
ggplot(merged,aes(-log10(V7),-log10(`p-value_gc`))) + geom_point()