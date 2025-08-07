version 1.2

# samtofastq
task samtofastq {
    input {
        File input_bam
        String prefix
        # Optional only if CRAM
        # File? reference_fasta
        # File? reference_fasta_index
        Int memoryMB
        Int num_cpu
        Int num_preempt
    }
    Int java_memory = floor(memoryMB - 500)
    Int diskGB = ceil(size(input_bam, "GB") * 7.5)

    String fastq1 = prefix + "_1.fastq"
    String fastq2 = prefix + "_2.fastq"
    String fastq0 = prefix + "_unpaired.fastq"
    String fastq1_gz = fastq1 + ".gz"
    String fastq2_gz = fastq2 + ".gz"

    String dollar = "$"

    command <<< 
        set -euxo pipefail

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Running SamToFastq..."

        java "-Xmx~{dollar}~{java_memory}m" -jar /opt/picard-tools/picard.jar SamToFastq \
            INPUT='~{input_bam}' \
            INCLUDE_NON_PF_READS=true \
            INCLUDE_NON_PRIMARY_ALIGNMENTS=false \
            VALIDATION_STRINGENCY=SILENT \
            FASTQ='~{fastq1}' \
            SECOND_END_FASTQ='~{fastq2}' \
            UNPAIRED_FASTQ='~{fastq0}'

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') SamToFastq completed successfully"

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Compressing FASTQs..."
        gzip -9 '~{fastq1}'
        gzip -9 '~{fastq2}'

        echo "Removing unpaired FASTQ"
        rm -f '~{fastq0}'

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Compression completed successfully"
    >>>

    output {
        File fastq1_out = fastq1_gz
        File fastq2_out = fastq2_gz
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
        preemptible: num_preempt
    }
}

# fastqc
task fastqc {
    input {
        File fastq1
        File fastq2
        Int memoryMB
        Int num_cpu
        Int num_preempt
    }

    Int diskGB = ceil((size(fastq1, "GB") + size(fastq2, "GB")) * 2)

    String fastq1_name = sub(sub(basename(fastq1), "\\.fastq.gz$", ""), "\\.fq.gz$", "")
    String fastq2_name = sub(sub(basename(fastq2), "\\.fastq.gz$", ""), "\\.fq.gz$", "")

    String fastq1_fastqc_html_str = fastq1_name + "_fastqc.html"
    String fastq1_fastqc_zip_str  = fastq1_name + "_fastqc.zip"
    String fastq1_fastqc_data_str = fastq1_name + ".fastqc_data.txt.gz"
    String fastq2_fastqc_html_str = fastq2_name + "_fastqc.html"
    String fastq2_fastqc_zip_str  = fastq2_name + "_fastqc.zip"
    String fastq2_fastqc_data_str = fastq2_name + ".fastqc_data.txt.gz"

    command <<< 
        set -euo pipefail
        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Running FastQC..."

        fastqc '~{fastq1}' '~{fastq2}' \
            --threads ~{num_cpu} \
            --outdir .

        unzip -p '~{fastq1_name}_fastqc.zip' '~{fastq1_name}_fastqc/fastqc_data.txt' | gzip > '~{fastq1_name}.fastqc_data.txt.gz'
        unzip -p '~{fastq2_name}_fastqc.zip' '~{fastq2_name}_fastqc/fastqc_data.txt' | gzip > '~{fastq2_name}.fastqc_data.txt.gz'

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') FastQC completed successfully"
    >>>

    output {
        File fastq1_fastqc_html = fastq1_fastqc_html_str
        File fastq1_fastqc_zip  = fastq1_fastqc_zip_str
        File fastq1_fastqc_data = fastq1_fastqc_data_str

        File fastq2_fastqc_html = fastq2_fastqc_html_str
        File fastq2_fastqc_zip  = fastq2_fastqc_zip_str
        File fastq2_fastqc_data = fastq2_fastqc_data_str
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
        preemptible: num_preempt
    }
}

# star
task star {

    input {
        File fastq1
        File fastq2
        String prefix
        File star_index

        # STAR args
        Int? outFilterMultimapNmax
        Int? alignSJoverhangMin
        Int? alignSJDBoverhangMin
        Int? outFilterMismatchNmax
        Float? outFilterMismatchNoverLmax
        Int? alignIntronMin
        Int? alignIntronMax
        Int? alignMatesGapMax
        String? outFilterType
        Float? outFilterScoreMinOverLread
        Float? outFilterMatchNminOverLread
        Int? limitSjdbInsertNsj
        String? outSAMstrandField
        String? outFilterIntronMotifs
        String? alignSoftClipAtReferenceEnds
        String? quantMode
        String? outSAMattrRGline
        String? outSAMattributes
        File? varVCFfile
        String? waspOutputMode
        Int? chimSegmentMin
        Int? chimJunctionOverhangMin
        String? chimOutType
        Int? chimMainSegmentMultNmax
        Int? chimOutJunctionFormat
        File? sjdbFileChrStartEnd
        String? quantTranscriptomeSAMoutput
        Int? winAnchorMultimapNmax
        String? genomeTransformOutput

        Int memoryMB
        Int num_cpu
        Int num_preempt
    }

    Int diskGB = ceil((size(fastq1, "GB") + size(fastq2, "GB") + size(star_index, "GB")) * 3)

    command <<<
        set -euo pipefail

        echo "FASTQs:"
        echo "~{fastq1}"
        echo "~{fastq2}"

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Extracting STAR index..."

        # Extract the star index
        mkdir -p star_index
        tar -xvvf "~{star_index}" -C star_index --strip-components=1

        # Place holders for optional outputs
        mkdir -p star_out
        touch star_out/"~{prefix}.Aligned.toTranscriptome.out.bam"
        touch star_out/"~{prefix}.Chimeric.out.sorted.bam"
        touch star_out/"~{prefix}.Chimeric.out.sorted.bam.bai"
        touch star_out/"~{prefix}.ReadsPerGene.out.tab"

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Running STAR alignment..."
        python3 /src/run_STAR.py \
            star_index '~{fastq1}' '~{fastq2}' '~{prefix}' \
            --output_dir star_out \
            ~{"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ~{"--alignSJoverhangMin " + alignSJoverhangMin} \
            ~{"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ~{"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ~{"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ~{"--alignIntronMin " + alignIntronMin} \
            ~{"--alignIntronMax " + alignIntronMax} \
            ~{"--alignMatesGapMax " + alignMatesGapMax} \
            ~{"--outFilterType " + outFilterType} \
            ~{"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ~{"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ~{"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ~{"--outSAMstrandField " + outSAMstrandField} \
            ~{"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ~{"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ~{"--quantMode " + quantMode} \
            ~{"--outSAMattrRGline " + outSAMattrRGline} \
            ~{"--outSAMattributes " + outSAMattributes} \
            ~{"--varVCFfile " + varVCFfile} \
            ~{"--waspOutputMode " + waspOutputMode} \
            ~{"--chimSegmentMin " + chimSegmentMin} \
            ~{"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ~{"--chimOutType " + chimOutType} \
            ~{"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ~{"--chimOutJunctionFormat " + chimOutJunctionFormat} \
            ~{"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            ~{"--quantTranscriptomeSAMoutput " + quantTranscriptomeSAMoutput} \
            ~{"--winAnchorMultimapNmax " + winAnchorMultimapNmax} \
            ~{"--genomeTransformOutput " + genomeTransformOutput} \
            --threads ~{num_cpu}

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') STAR alignment completed successfully"

    >>>

    output {
        File bam_file = "star_out/~{prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "star_out/~{prefix}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "star_out/~{prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "star_out/~{prefix}.Chimeric.out.junction.gz"
        File chimeric_bam_file = "star_out/~{prefix}.Chimeric.out.sorted.bam"
        File chimeric_bam_index = "star_out/~{prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "star_out/~{prefix}.ReadsPerGene.out.tab.gz"
        File junctions = "star_out/~{prefix}.SJ.out.tab.gz"
        File junctions_pass1 = "star_out/~{prefix}._STARpass1/~{prefix}.SJ.pass1.out.tab.gz"
        Array[File] logs = [
            "star_out/~{prefix}.Log.final.out",
            "star_out/~{prefix}.Log.out",
            "star_out/~{prefix}.Log.progress.out"
        ]
    }

    runtime {
        # V11 used 2.7.11a that doesn't have quantTranscriptomeSAMoutput
        # V11 uses 2.7.11b that does have quantTranscriptomeSAMoutput BanSingleEnd_ExtendSoftclip
        # https://github.com/alexdobin/STAR/releases
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11" 
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
        preemptible: num_preempt
    }
}

# markduplicates
task markduplicates {
    input {
        File input_bam
        String prefix
        Int? max_records_in_ram
        Float? sorting_collection_size_ratio

        Int memoryMB
        Int num_cpu
        Int num_preempt
    }
    
    Int java_memory = floor(memoryMB - 500)
    Int diskGB = ceil(size(input_bam, "GB") * 3)
    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")
    String output_bai = output_bam + ".bai"
    String metrics_file = prefix + ".marked_dup_metrics.txt"
    String dollar = "$"

    command <<<
        set -euo pipefail

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Running MarkDuplicates..."

        java "-Xmx~{dollar}~{java_memory}m" -jar /opt/picard-tools/picard.jar MarkDuplicates \
            I='~{input_bam}' \
            O='~{output_bam}' \
            M='~{metrics_file}' \
            PROGRAM_RECORD_ID=null \
            ASSUME_SORT_ORDER=coordinate \
            TMP_DIR=. \
            ~{"MAX_RECORDS_IN_RAM=" + max_records_in_ram} \
            ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio}

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') MarkDuplicates completed successfully"

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Indexing marked duplicate bam..."
        samtools index '~{output_bam}'
        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Indexing completed successfully"
    >>>

    output {
        File bam_file = output_bam
        File bam_index = output_bai
        File metrics = metrics_file
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
        preemptible: num_preempt
    }
}


# rsem
task rsem {
    input {
        File transcriptome_bam
        File rsem_reference
        String prefix
        Int? max_frag_len
        Boolean estimate_rspd = true
        Boolean is_stranded = false
        Boolean paired_end = true

        Int memoryMB
        Int num_cpu
        Int num_preempt
    }

    Int diskGB = ceil(size(transcriptome_bam, "GB") * 2.5)

    command <<< 
        set -euo pipefail

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Extracting RSEM reference..."
        mkdir rsem_ref
        tar -xvvf '~{rsem_reference}' -C rsem_ref --strip-components=1

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Running RSEM calculate-expression..."
        rsem-calculate-expression \
            --num-threads ~{num_cpu} \
            --no-bam-output \
            ~{"--fragment-length-max " + max_frag_len} \
            ~{if paired_end then "--paired-end" else ""} \
            ~{if estimate_rspd then "--estimate-rspd" else ""} \
            ~{if is_stranded then "--forward-prob 0.0" else ""} \
            --bam '~{transcriptome_bam}' \
            rsem_ref/rsem_reference \
            '~{prefix}.rsem'
        echo "$(date +'[%Y-%m-%d %H:%M:%S]') RSEM completed successfully"

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Compressing outputs..."
        gzip -f '~{prefix}.rsem.genes.results'
        gzip -f '~{prefix}.rsem.isoforms.results'
        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Compression completed successfully"
    >>>

    output {
        File genes = "~{prefix}.rsem.genes.results.gz"
        File isoforms = "~{prefix}.rsem.isoforms.results.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
        preemptible: num_preempt
    }
}

# rnaseqc2
task rnaseqc2 {
    input {
        File bam_file
        File genes_gtf
        String sample_id

        String? strandedness
        File? intervals_bed
        File? reference_fasta
        File? reference_fasta_index
        String? flags

        Int memoryMB
        Int num_cpu
        Int num_preempt
    }

    Int diskGB = ceil((size(bam_file, "GB") + size(genes_gtf, "GB") + size(intervals_bed, "GB") + size(reference_fasta, "GB")) * 2.5)

    command <<< 
        set -euo pipefail
        # Touch placeholder files
        touch '~{sample_id}.fragmentSizes.txt'
        touch '~{sample_id}.gc_content.tsv'

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Running RNA-SeQC 2..."
        rnaseqc '~{genes_gtf}' '~{bam_file}' . \
            -s '~{sample_id}' \
            ~{"--bed " + intervals_bed} \
            ~{"--stranded " + strandedness} \
            ~{"--fasta " + reference_fasta} \
            ~{"-vv " + flags}
        echo "$(date +'[%Y-%m-%d %H:%M:%S]') RNA-SeQC 2 completed successfully"

        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Compressing GCT outputs..."
        gzip -f *.gct
        echo "$(date +'[%Y-%m-%d %H:%M:%S]') Compression completed successfully"
        
    >>>

    output {
        File gene_tpm = "~{sample_id}.gene_tpm.gct.gz"
        File gene_counts = "~{sample_id}.gene_reads.gct.gz"
        File exon_counts = "~{sample_id}.exon_reads.gct.gz"
        File metrics = "~{sample_id}.metrics.tsv"
        File gc_content = "~{sample_id}.gc_content.tsv"
        File insertsize_distr = "~{sample_id}.fragmentSizes.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "~{memoryMB} MB"
        disks: "local-disk ~{diskGB} HDD"
        cpu: num_cpu
        preemptible: num_preempt
    }

}