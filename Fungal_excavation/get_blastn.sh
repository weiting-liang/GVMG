for i in `ls our_MAGs/*.fa`;do id=`basename $i .fa` && echo "blastn -query $i -db /path/database/NCBI/blast_nt/nt -evalue 1e-5 -outfmt 6 -num_threads 6 -out ./blastn/$id.out" >> blastn.sh;done
