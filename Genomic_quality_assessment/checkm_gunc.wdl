version 1.0

workflow batch_analysis {
    input {
        String Batch_id
        String input_dir
        String Gunc_database
        String docker_image = "stereonote_hpc_external/liangweiting_3f25d44e851f4c6e8edb6cad1d6ce377_private:latest"
    }

    call Unzip {
        input:
            input_dir = input_dir
    }

    call CheckM2 {
        input:
            Batch_id = Batch_id,
            input_dir = Unzip.unzipped_dir
    }

    call GUNC {
        input:
            Batch_id = Batch_id,
            input_dir = Unzip.unzipped_dir,
            gunc_db = Gunc_database,
            docker_image = docker_image
    }

    output {
        File checkm2_result = CheckM2.checkm2_result
        File gunc_result = GUNC.gunc_result
    }
}

task Unzip {
    input {
        String input_dir  # 保持input_dir为String类型
    }

    command {
        mkdir -p unzipped_dir
        mkdir -p temp_unzip_dir
        tar -xzvf ${input_dir} -C temp_unzip_dir
        mv temp_unzip_dir/*/* unzipped_dir
    }

    output {
        File unzipped_dir = "unzipped_dir"
    }

    runtime {
        docker_url: "stereonote_hpc_external/chenjunhong_c5a0b8b1b8514512a59605522102b530_private:latest"
        req_cpu: 1
        req_memory: "2Gi"
    }
}

task CheckM2 {
    input {
        String Batch_id
        File input_dir 
    }

    command {
        #checkm v1.0.2
        export PATH=/opt/conda/envs/checkm2/bin:$PATH
        export TMPDIR=/short/tmpdir 
        checkm2 predict -i ${input_dir} -o ${Batch_id}_checkm2_result --database_path /opt/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd -x fa -t 8 --force
    }

    output {
        File checkm2_result = "${Batch_id}_checkm2_result/quality_report.tsv"
    }

    runtime {
        docker_url: "stereonote_hpc_external/chenjunhong_c5a0b8b1b8514512a59605522102b530_private:latest"
        req_cpu: 8
        req_memory: "35Gi"
    }
} 

task GUNC {
    input {
        String Batch_id
        File input_dir
        String gunc_db
        String docker_image
    }

    command {
        export PATH=/opt/conda/envs/ProAssessment/bin:$PATH
        mkdir -p ${Batch_id}_gunc_result
        gunc run --db_file ${gunc_db} --input_dir ${input_dir} --file_suffix .fa --threads 8 --out_dir ${Batch_id}_gunc_result
    }

    output {
        File gunc_result = "${Batch_id}_gunc_result/GUNC.progenomes_2.1.maxCSS_level.tsv"
    }

    runtime {
        docker_url: docker_image
        req_cpu: 8
        req_memory: "60Gi"
    }
}
