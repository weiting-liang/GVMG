mkdir -p scaffolds.fa
for i in `cat scaffolds.fa.gz.list`;do id=`basename $i .gz` && gunzip -c $i > scaffolds.fa/$id ;done
for file in ./scaffolds.fa/*.fa; do
  mash sketch -o "${file%.fa}" $file
done

ls ./scaffolds.fa/*.msh > ./scaffolds.fa/query_list.txt
for i in `cat ./scaffolds.fa/query_list.txt `;do mash dist $i -l ./scaffolds.fa/query_list.txt >> mash_distances.tsv ;done
