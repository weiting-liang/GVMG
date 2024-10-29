version 1.0

workflow batch_analysis {
    input {
        String Batch_id
        String input_dir 
        String rRNA_5S_database
        String rRNA_16S_database
        String rRNA_23S_database
        String docker_image = "stereonote_hpc_external/liangweiting_3f25d44e851f4c6e8edb6cad1d6ce377_private:latest"
    }

    call Unzip {
        input:
            input_dir = input_dir
    }

    call Infernal {
        input:
            Batch_id = Batch_id,
            input_dir = Unzip.unzipped_dir,
            infernal_5S = rRNA_5S_database,
            infernal_16S = rRNA_16S_database,
            infernal_23S = rRNA_23S_database,
            docker_image = docker_image
    }

    output {
        File infernal_result = Infernal.infernal_result
    }
}

task Unzip {
    input {
        String input_dir
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

task Infernal {
    input {
        String Batch_id
        File input_dir
        String infernal_5S
        String infernal_16S
        String infernal_23S
        String docker_image
    }

    command {
        export PATH=/opt/conda/envs/ProAssessment/bin:$PATH
        for i in `ls ${input_dir}/*.fa`; do
            id=`basename $i .fa`
            cmsearch --tblout $id.5S_rna_results.tbl ${infernal_5S} $i
            cmsearch --tblout $id.16S_rna_results.tbl ${infernal_16S} $i
            cmsearch --tblout $id.23S_rna_results.tbl ${infernal_23S} $i
            result_5S="N"
            result_16S="N"
            result_23S="N"
            if grep -qv "^#" "$id.5S_rna_results.tbl"; then
                result_5S="Y"
            fi
            if grep -qv "^#" "$id.16S_rna_results.tbl"; then
                result_16S="Y"
            fi
            if grep -qv "^#" "$id.23S_rna_results.tbl"; then
                result_23S="Y"
            fi
            result_file="$id_5S_16S_23S.tsv"
            echo -e "$id\t$result_5S\t$result_16S\t$result_23S" > $id.5S_16S_23S.tsv
        done
        cat *.5S_16S_23S.tsv > ${Batch_id}.5S_16S_23S.rRNA.tsv
    }

    output {
        File infernal_result = "${Batch_id}.5S_16S_23S.rRNA.tsv"
    }

    runtime {
        docker_url: docker_image
        req_cpu: 2
        req_memory: "8Gi"
    }
}
