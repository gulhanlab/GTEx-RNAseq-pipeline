version development

workflow samtools_idx {
    input {
        File bam
    }

    call samtoolsIndex {
        input:
            bam = bam
    }

    output {
        File bai = samtoolsIndex.bai
    }
}


task samtoolsIndex {
    input {
        File bam
        Int memoryMB = 1024
        Int boot_diskGB = 8
        Int num_cpu = 1
        Int time = 30 
        Int num_preempt = 1
        Int max_retries = 1
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
        memory: memoryMB + " MB"
        cpu: num_cpu
        runtime_minutes: time
        bootDiskSizeGb: boot_diskGB
        disks: "local-disk " + diskGB + " HDD"
        preemptible: num_preempt
        maxRetries: max_retries
    }
}