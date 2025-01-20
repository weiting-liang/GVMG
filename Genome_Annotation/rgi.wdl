version 1.0

workflow run_rgi {
    input {
        Array[String] fa_files
    }

    call RunRgi {
        input:
            fa_files = fa_files
    }

    output {
        File compressed_results = RunRgi.compressed_results
    }
}

task RunRgi {
    input {
        Array[String] fa_files
    }

    command {
        export PATH=/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv3.7/bin:$PATH
        export PATH=/hwfssz1/ST_HEALTH/P18Z10200N0127/juyanmei/bin/miniconda3/envs/rgi/bin:$PATH
        mkdir -p rgi_result
        cp -r /hwfssz1/ST_HEALTH/P18Z10200N0127/juyanmei/database_tmp/card/* .
        echo "${sep='\n' fa_files}" > input.list
        while IFS= read -r fa_file; do
            prodigal -i $fa_file -a $(basename $fa_file).faa -p meta
            rgi main --input_sequence $(basename $fa_file).faa --output_file rgi_result/$(basename $fa_file) --local --clean -t protein -n 4
        done < input.list
        tar -czvf rgi_results.tar.gz rgi_result
    }

    output {
        File compressed_results = "rgi_results.tar.gz"
    }

    runtime {
        req_cpu: 4
        req_memory: "15Gi"
    }
}
