version development

import "tasks.wdl" as tasks

workflow rnaseq_pipeline_bam_workflow {

    input {
        String prefix
        File input_bam
        # File? reference_fasta
        # File? reference_fasta_index

        # SAMTOFASTQ
        Int samtofastq_memoryMB
        Int samtofastq_num_cpu
        Int samtofastq_num_preempt

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
        File genes_gtf
        File? intervals_bed
        Int rnaseqc2_memoryMB
        Int rnaseqc2_num_cpu
        Int rnaseqc2_num_preempt
    }

    call tasks.samtofastq {
        input:
            input_bam = input_bam,
            prefix = prefix,
            memoryMB = samtofastq_memoryMB,
            num_cpu = samtofastq_num_cpu,
            num_preempt = samtofastq_num_preempt
    }

    call tasks.star {
        input:
            fastq1 = samtofastq.fastq1_out,
            fastq2 = samtofastq.fastq2_out,
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
            genes_gtf = genes_gtf,
            sample_id = prefix,
            # reference_fasta = reference_fasta,
            # reference_fasta_index = reference_fasta_index,
            intervals_bed = intervals_bed,
            memoryMB = rnaseqc2_memoryMB,
            num_cpu = rnaseqc2_num_cpu,
            num_preempt = rnaseqc2_num_preempt
    }
    
    output {
        File samtofastq_fastq1 = samtofastq.fastq1_out
        File samtofastq_fastq2 = samtofastq.fastq2_out

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