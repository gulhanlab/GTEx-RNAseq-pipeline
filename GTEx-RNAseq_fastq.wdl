version development

import "tasks.wdl" as tasks

workflow rnaseq_pipeline_fastq_workflow {

    input {
        String prefix
        File input_fastq1
        File input_fastq2
        # File? reference_fasta
        # File? reference_fasta_index

        # FASTQC
        Int fastqc_memoryMB
        Int fastqc_num_cpu
        Int fastqc_num_preempt

        # STAR
        File star_index
        Int star_memoryMB
        Int star_num_cpu
        Int star_num_preempt

        # MARKDUPLICATES
        Int markduplicates_memoryMB
        Int markduplicates_num_cpu
        Int markduplicates_num_preempt

        # RSEM
        File rsem_reference
        Int rsem_memoryMB
        Int rsem_num_cpu
        Int rsem_num_preempt

        # RNASEQC2
        File rnaseqc2_genes_gtf
        File? rnaseqc2_intervals_bed
        Int rnaseqc2_memoryMB
        Int rnaseqc2_num_cpu
        Int rnaseqc2_num_preempt
    }

    call tasks.fastqc {
        input:
            fastq1 = input_fastq1,
            fastq2 = input_fastq2,
            memoryMB = fastqc_memoryMB,
            num_cpu = fastqc_num_cpu,
            num_preempt = fastqc_num_preempt
    }

    call tasks.star {
        input:
            fastq1 = input_fastq1,
            fastq2 = input_fastq2,
            prefix = prefix,
            star_index = star_index,
            memoryMB = star_memoryMB,
            num_cpu = star_num_cpu,
            num_preempt = star_num_preempt
    }

    call tasks.markduplicates {
        input:
            input_bam = star.bam_file,
            prefix = prefix,
            memoryMB = markduplicates_memoryMB,
            num_cpu = markduplicates_num_cpu,
            num_preempt = markduplicates_num_preempt
    }

    call tasks.rsem {
        input:
            transcriptome_bam = star.transcriptome_bam,
            rsem_reference = rsem_reference,
            prefix = prefix,
            memoryMB = rsem_memoryMB,
            num_cpu = rsem_num_cpu,
            num_preempt = rsem_num_preempt
    }

    call tasks.rnaseqc2 {
        input:
            bam_file = markduplicates.bam_file,
            genes_gtf = rnaseqc2_genes_gtf,
            sample_id = prefix,
            # reference_fasta = reference_fasta,
            # reference_fasta_index = reference_fasta_index,
            intervals_bed = rnaseqc2_intervals_bed,
            memoryMB = rnaseqc2_memoryMB,
            num_cpu = rnaseqc2_num_cpu,
            num_preempt = rnaseqc2_num_preempt
    }
    
    output {
        File fastqc1_html = fastqc.fastq1_fastqc_html
        File fastqc1_zip = fastqc.fastq1_fastqc_zip
        File fastqc1_data = fastqc.fastq1_fastqc_data
        File fastqc2_html = fastqc.fastq2_fastqc_html
        File fastqc2_zip = fastqc.fastq2_fastqc_zip
        File fastqc2_data = fastqc.fastq2_fastqc_data

        File star_bam = star.bam_file
        File star_bam_index = star.bam_index
        File star_transcriptome_bam = star.transcriptome_bam
        File star_chimeric_junctions = star.chimeric_junctions
        File star_chimeric_bam = star.chimeric_bam_file
        File star_chimeric_index = star.chimeric_bam_index
        File star_read_counts = star.read_counts
        File star_junctions = star.junctions
        File star_junctions_pass = star.junctions_pass1
        Array[File] star_logs = star.logs

        File markdup_bam = markduplicates.bam_file
        File markdup_bam_index = markduplicates.bam_index
        File markdup_metrics = markduplicates.metrics

        File rsem_genes = rsem.genes
        File rsem_isoforms = rsem.isoforms

        File rnaseqc2_gene_tpm = rnaseqc2.gene_tpm
        File rnaseqc2_gene_counts = rnaseqc2.gene_counts
        File rnaseqc2_exon_counts = rnaseqc2.exon_counts
        File rnaseqc2_metrics = rnaseqc2.metrics
        File rnaseqc2_gc_content = rnaseqc2.gc_content
        File rnaseqc2_insertsize = rnaseqc2.insertsize_distr
  }
}