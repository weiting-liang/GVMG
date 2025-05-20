#haolilan
# Let's continue to generate bracken distribution for profile correction
kraken2 --db=./ --threads=4 <( find -L library \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) -exec cat {} + ) >   database.kraken

#conda activate bracken
for l in 50 100 200 300; do
kmer2read_distr \
        --seqid2taxid seqid2taxid.map \
        --taxonomy ./taxonomy \
        --kraken database.kraken \
        --output "database${l}mers.kraken" \
        -l $l -t 1
    
    generate_kmer_distribution.py \
        -i "database${l}mers.kraken" \
        -o "database${l}mers.kmer_distrib"
done
 
# rm -r ../kraken_library/
