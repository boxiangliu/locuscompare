for i in {1..22} X Y; do 
bcftools query -f \
"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/EAS_AF\t%INFO/AMR_AF\t%INFO/AFR_AF\t%INFO/EUR_AF\t%INFO/SAS_AF\n" \
-r $i \
-H \
/mnt/data/shared/1KG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz \
> /srv/persistent/bliu2/locuscompare/data/1KG/GRCh37/chr$i.txt &
done 
