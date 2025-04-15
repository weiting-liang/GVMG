mkdir -p output/bin
mkdir -p report/multisample20
mkdir -p tmp/metaspades
mkdir -p tmp/metabat2

rm bin_Metabat2_se.sh

while read line
do
        id=`basename $line .rmhost.fq.gz`

        ##binning prepare
        echo "mkdir -p ./tmp/metabat2/$id" >> bin_Metabat2_se.sh
        echo "bwa index ./output/assembly/$id.contigs.fa.gz -p ./tmp/metabat2/$id/$id.contigs.fa.gz 2> ./log/$id.contigs.fa.gz.index.log ##make index" >> bin_Metabat2_se.sh
        echo "for i in \`cat ./multi20_list/$id.multi20.bwa.txt\`;do bwa mem -t 4 ./tmp/metabat2/$id/$id.contigs.fa.gz ./output/rmhost/$id.rmhost.fq.gz  2>> ./log/$id.contigs.fa.gz.align.log|samtools sort -m 3G -@4 -T ./tmp/metabat2/$id/\$i.align2_$id.contigs.sorted.bam -O BAM -o ./tmp/metabat2/$id/\$i.align2_$id.contigs.sorted.bam - 2>> ./log/$id.contigs.fa.gz.align.log && samtools index ./tmp/metabat2/$id/\$i.align2_$id.contigs.sorted.bam && samtools flagstat ./tmp/metabat2/$id/\$i.align2_$id.contigs.sorted.bam > ./tmp/metabat2/$id/\$i.align2_$id.contigs.sorted.bam.flagstat;done" >> bin_Metabat2_se.sh
        echo "bam_list=\`ls ./tmp/metabat2/$id/*.align2*.contigs.sorted.bam | xargs\` && jgi_summarize_bam_contig_depths --outputGC ./tmp/metabat2/$id/$id.outputGC --outputDepth ./tmp/metabat2/$id/$id.contigs.depth.txt \$bam_list 1> ./log/$id.coverage.log 2> ./log/$id.coverage.error" >> bin_Metabat2_se.sh

        ##binning
        echo "mkdir -p ./output/bin/$id" >> bin_Metabat2_se.sh
        echo "metabat2 -t 4 --minContig 1500 -i ./output/assembly/$id.contigs.fa.gz -a ./tmp/metabat2/$id/$id.contigs.depth.txt -o ./output/bin/$id/$id -v 2> ./log/$id.metabat2.log" >> bin_Metabat2_se.sh
        echo "rm -r ./tmp/metabat2/$id/" >> bin_Metabat2_se.sh

        echo "###########################################   sample $id finished !!! ################################################################################################" >> bin_Metabat2_se.sh
done < ./samples_rmhost1.txt
