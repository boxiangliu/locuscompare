gtex_fn=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_FOR_QC_ONLY
out_fn=/srv/persistent/bliu2/locuscompare/processed_data/count_gtex_samples/count.txt

for f in $(ls $gtex_fn/*.bed); do
    echo $f
    n=$(zcat $f | awk '{print NF}' | head -n1)
    n=$((n-4))
    echo $n
    echo -e "$f\t$n" >> $out_fn
done
