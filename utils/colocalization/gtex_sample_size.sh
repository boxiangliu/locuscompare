expr_dir=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/expression_matrices
out_fn=/srv/persistent/bliu2/locuscompare/processed_data/colocalization/gtex_sample_size/sample_size.txt
mkdir -p $(dirname $out_fn)

echo -e 'tissue\tsize' > $out_fn
for f in $(ls $expr_dir/*.gz); do
	tissue=$(basename $f | sed 's/.v7.normalized_expression.bed.gz//')
	n=$(zcat $f | head -n1 | awk '{print NF}')
	n=$((n-3))
	echo -e $tissue'\t'$n >> $out_fn
done