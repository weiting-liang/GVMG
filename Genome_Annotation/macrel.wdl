version 1.0

workflow run_macrel {
    input {
        Array[String] fa_files
    }

    call RunMacrel {
        input:
            fa_files = fa_files
    }

    output {
        File compressed_results = RunMacrel.compressed_results
    }
}

task RunMacrel {
    input {
        Array[String] fa_files
    }

    command {
        export PATH=/opt/software/miniconda3/envs/macrel_env/bin:$PATH
        mkdir -p ./tmp
        echo "${sep='\n' fa_files}" > input.list

        while IFS= read -r fa_file; do
            macrel contigs --fasta $fa_file --output ./output -t 4 --tag $(basename $fa_file .fa) --tmpdir ./tmp --force --log-file log.txt --log-append
        done < input.list

        mkdir AMP_predict
        mv output/*percontigs* AMP_predict/
        mv output/*prediction* AMP_predict/
        tar -czvf AMP_predict.tar.gz AMP_predict/
    }

    output {
        File compressed_results = "AMP_predict.tar.gz"
    }

    runtime {
        req_cpu:4
        req_memory: "20Gi"
        docker_url: "stereonote_hpc/liangweiting_a289de442c414982ae521309dec79d1e_private:latest"
    }
}
