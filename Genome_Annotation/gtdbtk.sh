#classify
export PATH=/path/anaconda3/envs/gtdbtk-2.4.0/bin:$PATH && export GTDBTK_DATA_PATH=/path/anaconda3/envs/gtdbtk-2.4.0/share/gtdbtk-2.4.0/db  && gtdbtk classify_wf --genome_dir /path/SGB --out_dir ./output --extension fa --skip_ani_screen --cpus 8 --pplacer_cpus 1 > ./log.txt

#phylogenetic tree
gtdbtk infer --msa_file gtdbtk.bac120.user_msa.fasta --out_dir PhyloTree --cpus 8 --prefix bac
