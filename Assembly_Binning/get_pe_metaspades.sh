#!/bin/bash
mkdir -p log
mkdir -p output/rmhost
mkdir -p output/assembly
mkdir -p tmp/metaspades

rm assembly_pe.sh

while read line
do
        id=`basename $line .rmhost.1.fq.gz`
        dir=`dirname $line`
        ##rmhost
        echo "bowtie2 --end-to-end --very-sensitive -p 8 -x /ldfssz1/ST_META/share/User/zhujie/database/galaxy_indexes/chm13v2/bowtie2_index/chm13v2 -1 $dir/$id.rmhost.1.fq.gz -2 $dir/$id.rmhost.2.fq.gz 2> ./log/$id.bowtie2.log | samtools fastq -N -c 5 -f 12 -1 ./output/rmhost/$id.rmhost.1.fq.gz -2 ./output/rmhost/$id.rmhost.2.fq.gz" >> assembly_pe.sh

        ##assembly
        echo "metaspades.py -1 ./output/rmhost/$id.rmhost.1.fq.gz -2 ./output/rmhost/$id.rmhost.2.fq.gz -k '21', '33', '55', '77' --memory 81 --threads 8 --checkpoints last -o ./tmp/metaspades/$id > ./log/$id.metaspades.log && pigz -p 8 ./tmp/metaspades/$id/scaffolds.fasta && mv ./tmp/metaspades/$id/scaffolds.fasta.gz ./output/assembly/$id.scaffolds.fasta.gz && pigz -p 8 ./tmp/metaspades/$id/contigs.fasta && mv ./tmp/metaspades/$id/contigs.fasta.gz ./output/assembly/$id.contigs.fasta.gz && rm -r ./tmp/metaspades/$id/" >> assembly_pe.sh

done < ./samples_rmhost1.txt
