mkdir -p log
mkdir -p output/rmhost
mkdir -p output/assembly
mkdir -p tmp/megahit

rm assembly_se.sh

while read line
do
        id=`basename $line .rmhost.fq.gz`
        dir=`dirname $line`
        ##rmhost
        echo "bowtie2 --end-to-end --very-sensitive -p 4 -x /ldfssz1/ST_META/share/User/zhujie/database/galaxy_indexes/chm13v2/bowtie2_index/chm13v2 -U $dir/$id.rmhost.fq.gz 2> ./log/$id.bowtie2.log | samtools fastq -N -c 5 -f 4 -F 256 --threads 4 | pigz -c -p 4 > ./output/rmhost/$id.rmhost.fq.gz" >> assembly_se.sh
        ##assembly
        echo "megahit -r ./output/rmhost/$id.rmhost.fq.gz -t 4 --presets meta-sensitive --min-contig-len 100 --out-dir ./tmp/megahit/$id --out-prefix $id 2> ./log/$id.megahit.log && pigz -p 4 ./tmp/megahit/$id/$id.contigs.fa && mv ./tmp/megahit/$id/$id.contigs.fa ./output/assembly/$id.contigs.fa.gz && rm -r ./tmp/megahit/$id/" >> assembly_se.sh
done < ./samples_rmhost1.txt
