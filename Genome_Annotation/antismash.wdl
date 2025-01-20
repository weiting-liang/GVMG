version 1.0

workflow batch_analysis {
    input {
        Array[String] genome_ids
        Array[File] genomes
    }

    scatter (i in range(length(genome_ids))) {
        call RunAntismash {
            input:
                genome_id = genome_ids[i],
                genome = genomes[i]
        }
    }

    output {
        Array[File] compressed_results = RunAntismash.compressed_result
    }
}

task RunAntismash {
    input {
        String genome_id
        File genome
    }

    command {
        export PATH=/opt/software/miniconda3/envs/antismash_env/bin:$PATH
        antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --output-dir ./bgc_${genome_id} ${genome} --genefinding-tool prodigal -c 4
        tar -czvf ${genome_id}_bgc_result.tar.gz bgc_${genome_id}
    }

    output {
        File compressed_result = "${genome_id}_bgc_result.tar.gz"
    }

    runtime {
        docker_url: "stereonote_hpc/liangweiting_347567e0bc544cb4852d6b74d16cd1ea_private:latest"
        req_cpu: 4
        req_memory: "20Gi"
    }
}
