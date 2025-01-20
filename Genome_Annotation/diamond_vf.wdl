version 1.0

workflow batch_analysis {
    input {
        Array[String] genome_ids
        Array[File] genomes
        String vfdb
    }

    scatter (i in range(length(genome_ids))) {
        call RunFunctionalAnalysis {
            input:
                genome_id = genome_ids[i],
                genome = genomes[i],
                vfdb = vfdb
        }
    }

    output {
        Array[File] compressed_results = RunFunctionalAnalysis.compressed_result
    }
}

task RunFunctionalAnalysis {
    input {
        String genome_id
        File genome
        String vfdb
    }

    command {
        export PATH=/opt/software/miniconda3/envs/func_env/bin:$PATH
        prodigal -i ${genome} -a ${genome_id}.faa -p meta
        diamond blastp --db ${vfdb} -q ${genome_id}.faa -o vf_${genome_id}.matches_fmt6.txt --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
        tar -czvf ${genome_id}_functional2_result.tar.gz ${genome_id}.faa vf_${genome_id}.matches_fmt6.txt
    }

    output {
        File compressed_result = "${genome_id}_functional2_result.tar.gz"
    }

    runtime {
        docker_url: "stereonote_hpc/liangweiting_da8b4a29276a474a909f8f18e2359dc7_private:latest"
        req_cpu: 2
        req_memory: "5Gi"
    }
}
