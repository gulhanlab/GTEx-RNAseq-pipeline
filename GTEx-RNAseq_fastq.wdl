# WDL draft-2

import "tasks.wdl" as tasks

workflow rnaseq_pipeline_fastq_workflow {
    String prefix

    call tasks.fastqc {}

    call tasks.star {
        input: 
            prefix = prefix
    }

    call tasks.markduplicates {
        input: 
            input_bam = star.bam_file, 
            prefix = prefix
    }

    call tasks.rsem {
        input: 
            transcriptome_bam = star.transcriptome_bam, 
            prefix = prefix
    }

    call tasks.rnaseqc2 {
        input: 
            bam_file = markduplicates.bam_file, 
            sample_id = prefix
    }
}
