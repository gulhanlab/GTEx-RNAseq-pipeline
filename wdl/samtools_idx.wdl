version development

workflow samtools_index_workflow {
    input {
        File bam
    }

    call samtools_index {
        input:
            bam = bam
    }

    output {
        File bai = samtools_index.bai
    }
}


task samtools_index {
    input {
        File bam
        Int memoryMB = 1024
        Int boot_diskGB = 8
        Int num_cpu = 1
        String docker_image = "staphb/samtools:1.21"
    }
    # COMPUTE DISK SIZE
    Int diskGB = ceil(1.2 * size(bam, "GB"))
    
    command <<<
        samtools index ~{bam}
    >>>

    output {
        # By default, samtools creates an index file named <bam>.bai.
        File bai = "~{bam}.bai"
    }
    
    runtime {
        docker: docker_image
        bootDiskSizeGb: boot_diskGB
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
    }
}