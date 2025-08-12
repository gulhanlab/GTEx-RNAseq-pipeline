version development

workflow star_index_workflow {
    input {
        String star_index_name
        File reference_fasta
        File reference_fasta_index
        File gencode_annotation_gtf
        Int overhang # fastq read length - 1
        Boolean? use_hg19

        Int boot_disk_sizeGB = 10
        Int memoryMB = 49152
        Int num_cpu = 8
    }

    String docker_image = if defined(use_hg19) then "gulhanlab/gtex-rnaseq-pipeline:hg19_v2" else "gulhanlab/gtex-rnaseq-pipeline:hg38_v2"

    call star_index {
        input:
            star_index_name = star_index_name,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            gencode_annotation_gtf = gencode_annotation_gtf,
            overhang = overhang,
            boot_disk_sizeGB = boot_disk_sizeGB,
            memoryMB = memoryMB,
            num_cpu = num_cpu,
            docker_image = docker_image
    }

    output {
        File output_star_index = star_index.star_index_tar
    }
}


task star_index {
    input {
        String star_index_name
        File reference_fasta
        File reference_fasta_index
        File gencode_annotation_gtf
        Int overhang # fastq read length - 1

        Int boot_disk_sizeGB
        Int memoryMB = 49152
        Int num_cpu = 8
        
        String? docker_image = "gulhanlab/gtex-rnaseq-pipeline:hg38_v2" 
    }

    Int diskGB = 150
    
    command <<<
        set -euo pipefail

        # Create STAR index directory
        mkdir -p ~{star_index_name}
        STAR \
            --runMode genomeGenerate \
            --genomeDir ~{star_index_name} \
            --genomeFastaFiles ~{reference_fasta} \
            --sjdbGTFfile ~{gencode_annotation_gtf} \
            --sjdbOverhang ~{overhang} \
            --runThreadN ~{num_cpu}

        # Create a tarball of the STAR index
        tar -cvzf ~{star_index_name}.tar.gz ~{star_index_name}
    >>>

    output {
        File star_index_tar = "~{star_index_name}.tar.gz"
    }
    
    runtime {
        docker: docker_image
        bootDiskSizeGb: boot_disk_sizeGB
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
    }
}