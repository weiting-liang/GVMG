mkdir -p output/bin
mkdir -p report/multisample20
mkdir -p tmp/metaspades
mkdir -p tmp/metabat2

rm bin_Metabat2_pe.sh

while read line
do
        id=`basename $line .rmhost.1.fq.gz`

        ##binning prepare
        echo "mkdir -p ./tmp/metabat2/$id" >> bin_Metabat2_pe.sh
        echo "bwa index ./output/assembly/$id.scaffolds.fasta.gz -p ./tmp/metabat2/$id/$id.scaffolds.fasta.gz 2> ./log/$id.scaffolds.fasta.gz.index.log ##make index" >> bin_Metabat2_pe.sh
        echo "for i in \`cat ./multi20_list/$id.multi20.bwa.txt\`;do bwa mem -t 4 ./tmp/metabat2/$id/$id.scaffolds.fasta.gz -1 ./output/rmhost/$id.rmhost.1.fq.gz -2 ./output/rmhost/$id.rmhost.2.fq.gz  2>> ./log/$id.scaffolds.fasta.gz.align.log|samtools sort -m 3G -@4 -T ./tmp/metabat2/$id/\$i.align2_$id.scaffolds.sorted.bam -O BAM -o ./tmp/metabat2/$id/\$i.align2_$id.scaffolds.sorted.bam - 2>> ./log/$id.scaffolds.fa.gz.align.log && samtools index ./tmp/metabat2/$id/\$i.align2_$id.scaffolds.sorted.bam && samtools flagstat ./tmp/metabat2/$id/\$i.align2_$id.scaffolds.sorted.bam > ./tmp/metabat2/$id/\$i.align2_$id.scaffolds.sorted.bam.flagstat;done" >> bin_Metabat2_spe.sh
        echo "bam_list=\`ls ./tmp/metabat2/$id/*.align2*.scaffolds.sorted.bam | xargs\` && jgi_summarize_bam_contig_depths --outputGC ./tmp/metabat2/$id/$id.outputGC --outputDepth ./tmp/metabat2/$id/$id.scaffolds.depth.txt \$bam_list 1> ./log/$id.scaffolds.log 2> ./log/$id.coverage.error" >> bin_Metabat2_pe.sh

        ##binning
        echo "mkdir -p ./output/bin/$id" >> bin_Metabat2_pe.sh
        echo "metabat2 -t 4 --minContig 1500 -i ./output/assembly/$id.scaffolds.fa.gz -a ./tmp/metabat2/$id/$id.scaffolds.depth.txt -o ./output/bin/$id/$id -v 2> ./log/$id.metabat2.log" >> bin_Metabat2_pe.sh
        echo "rm -r ./tmp/metabat2/$id/" >> bin_Metabat2_pe.sh

        echo "###########################################   sample $id finished !!! ################################################################################################" >> bin_Metabat2_pe.sh
done < ./samples_rmhost1.txt
