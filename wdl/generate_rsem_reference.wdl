version development

workflow rsem_reference_workflow {
    input {
        String rsem_reference_name
        File reference_fasta
        File reference_fasta_index
        File gencode_annotation_gtf
        Boolean? use_hg19

        Int boot_disk_sizeGB = 10
        Int memoryMB = 4096
        Int num_cpu = 1
    }

    String docker_image = if defined(use_hg19) then "gulhanlab/gtex-rnaseq-pipeline:hg19_v2" else "gulhanlab/gtex-rnaseq-pipeline:hg38_v2"

    call rsem_reference {
        input:
            rsem_reference_name = rsem_reference_name,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            gencode_annotation_gtf = gencode_annotation_gtf,
            boot_disk_sizeGB = boot_disk_sizeGB,
            memoryMB = memoryMB,
            num_cpu = num_cpu,
            docker_image = docker_image
    }

    output {
        File output_rsem_reference = rsem_reference.rsem_reference_tar
    }
}


task rsem_reference {
    input {
        String rsem_reference_name
        File reference_fasta
        File reference_fasta_index
        File gencode_annotation_gtf

        Int boot_disk_sizeGB
        Int memoryMB = 4096
        Int num_cpu = 1
        
        String? docker_image = "gulhanlab/gtex-rnaseq-pipeline:hg38_v2" 
    }

    Int diskGB = 50
    
    command <<<
        set -euo pipefail

        mkdir -p ~{rsem_reference_name}
        cd ~{rsem_reference_name}

        # Prepare RSEM reference folder named rsem_reference_name
        # Files within this folder will have a prefix of "rsem_reference"
        rsem-prepare-reference \
            ~{reference_fasta} \
            rsem_reference \
            --gtf ~{gencode_annotation_gtf} \
            --num-threads ~{num_cpu}

        cd ..
        tar -cvzf ~{rsem_reference_name}.tar.gz ~{rsem_reference_name}
    >>>

    output {
        # Extracting the tar will contain the following files:
        # ~{rsem_reference_name}/rsem_reference.chrlist
        # ~{rsem_reference_name}/rsem_reference.grip
        # ~{rsem_reference_name}/rsem_reference.idx.fa
        # ~{rsem_reference_name}/rsem_reference.n2g.idx.fa
        # ~{rsem_reference_name}/rsem_reference.seq
        # ~{rsem_reference_name}/rsem_reference.ti
        # ~{rsem_reference_name}/rsem_reference.transcripts.fa
        File rsem_reference_tar = "~{rsem_reference_name}.tar.gz"
    }

    runtime {
        docker: docker_image
        bootDiskSizeGb: boot_disk_sizeGB
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
    }
}