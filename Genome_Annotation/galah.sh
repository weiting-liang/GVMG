#galah 0.4.0
/hwfssz1/ST_HEALTH/P18Z10200N0127/juyanmei/bin/miniconda3/bin/galah cluster \
  --genome-fasta-list Pro_mags_genomes.all.path.list \
  --checkm2-quality-report checkm2_report.txt \
  --genome-fasta-extension fa \
  --ani 95 \
  --precluster-ani 90 \
  --precluster-method finch \
  --min-completeness 50 \
  --max-contamination 5 \
  --fragment-length 1500 \
  --min-aligned-fraction 50 \
  --output-representative-fasta-directory ./hmq_genome_rep_cluster_species \
  --output-cluster-definition ./hmq_genome_rep_cluster_definition_species.tsv \
  --output-representative-list ./hmq_genome_rep_cluster_list_species.tsv \
  --threads 16 > hmq_genome_galah_cluster_species.log 2>&1
