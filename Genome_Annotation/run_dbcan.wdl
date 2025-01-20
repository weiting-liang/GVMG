version 1.0

workflow run_dbcan {
    input {
        Array[String] fa_files
    }

    call RunDbcan {
        input:
            fa_files = fa_files
    }

    output {
        File compressed_results = RunDbcan.compressed_results
    }
}

task RunDbcan {
    input {
        Array[String] fa_files
    }

    command {
        export PATH=/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv3.7/bin:$PATH
        export PATH=/ldfssz1/ST_META/share/User/zhujie/.conda/envs/dbcan/bin:$PATH
        mkdir -p dbcan_result
        echo "${sep='\n' fa_files}" > input.list
        while IFS= read -r fa_file; do
            prodigal -i $fa_file -a $(basename $fa_file).faa -p meta
            run_dbcan.py $(basename $fa_file).faa protein --db_dir /hwfssz1/ST_HEALTH/P18Z10200N0127/liangweiting/program/project_vagina2/Assay/06_function/cazy_dbcan/db --out_dir dbcan_result/$(basename $fa_file)
        done < input.list
        tar -czvf dbcan_results.tar.gz dbcan_result
    }

    output {
        File compressed_results = "dbcan_results.tar.gz"
    }

    runtime {
        req_cpu:2
        req_memory: "15Gi"
        docker_url: "stereonote_hpc/liangweiting_da8b4a29276a474a909f8f18e2359dc7_private:latest"
    }
}
