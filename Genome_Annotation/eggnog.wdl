version 1.0

workflow batch_analysis {
    input {
        Array[String] genome_ids
        Array[File] genomes
    }

    scatter (i in range(length(genome_ids))) {
        call RunEggnog {
            input:
                genome_id = genome_ids[i],
                genome = genomes[i]
        }
    }

    output {
        Array[File] compressed_results = RunEggnog.compressed_result
    }
}

task RunEggnog {
    input {
        String genome_id
        File genome
    }

    command {
        export PATH=/ldfssz1/ST_META/share/User/juyanmei/miniconda3/bin:$PATH
        mkdir result tmp
        prodigal -i ${genome} -a ${genome_id}.faa -p meta
        /ldfssz1/ST_META/share/User/juyanmei/miniconda3/bin/emapper.py --cpu 8 --data_dir /hwfssz1/ST_HEALTH/P18Z10200N0127/juyanmei/database_tmp/emapperdb_tmp/emapperdb-5.0.2 -o ${genome_id} --output_dir ./result --temp_dir ./tmp --override -i ${genome_id}.faa --dmnd_db /hwfssz1/ST_HEALTH/P18Z10200N0127/juyanmei/database_tmp/emapperdb_tmp/emapperdb-5.0.2/eggnog_proteins.dmnd --no_file_comments
        tar -czvf kegg_${genome_id}.tar.gz result
    }

    output {
        File compressed_result = "kegg_${genome_id}.tar.gz"
    }

    runtime {
        docker_url: "stereonote_hpc/liangweiting_da8b4a29276a474a909f8f18e2359dc7_private:latest"
        req_cpu: 8
        req_memory: "30Gi"
    }
}
