ls ./mags_fa_gz/*.fa.gz > list.txt

`rm work.sh`
while read line
do
    name=`basename $line .fa.gz`
    echo "###one sample start#####
perl ./getlength.pl $line $name.len
/path/blast/filter_blast  -i ../blastn/$name.out  -o $name.out.filter  --identity  90 --qfile $name.len --qper 80 --tops 20
grep '>' ../our_MAGs/$name.fa | wc -l > $name.contigsCount
grep '>' ../our_MAGs/$name.fa |cut -d ' ' -f 1 | sed 's/>//g' | xargs | tr ' ' ',' > $name.contigs
echo '$name.fa'> $name.txt
paste $name.txt $name.contigsCount $name.contigs > $name.cluster
cut -f 2 $name.out.filter > $name.target.list

for i in \`cat $name.target.list\`;do grep \$i /path/accession2taxid/nucl_gb.accession2taxid >> $name.target_taxid.txt;done
sort $name.target_taxid.txt |uniq > $name.target_taxid2.txt
cp $name.out.filter $name.out.filter.rename
cat $name.target_taxid2.txt |while read line;do target=\`echo \$line | cut -d \" \" -f 2\` && taxid=\`echo \$line | cut -d \" \" -f 3\` && echo \"sed -i 's/\$target/\$taxid/g' $name.out.filter.rename\" >> $name.rename2.sh;done
sh $name.rename2.sh
awk -v OFS='\t' '{print \$1,\$3,\$2}' $name.out.filter.rename > $name.tax
perl ./mlg.tax.perl $name.tax $name.cluster > $name.taxo
####one sample ending##############
" >> work.sh

done < list.txt
