# WDL draft-2

# samtofastq
task samtofastq {

    File input_bam
    String prefix
    # File? reference_fasta
    # File? reference_fasta_index

    Float memory
    Int java_memory = floor(memory - 0.5)
    # Int disk_space
    Int num_threads
    Int num_preempt

    Int diskGB = ceil(size(input_bam) * 5) 

    command {
        set -euo pipefail

        # make sure path is absolute
        input_bam_abs=${input_bam}
        if [[ $input_bam_abs != /* ]]; then
            input_bam_abs=$PWD/$input_bam_abs
        fi

        mkdir samtofastq  # workaround for named pipes
        echo $(date +"[%Y-%m-%d %H:%M:%S] Running SamToFastq...")
        python3 -u /src/run_SamToFastq.py $input_bam_abs -p ${prefix} --output_dir samtofastq --memory ${java_memory}
        mv samtofastq/${prefix}_*.fastq.gz .
        echo $(date +"[%Y-%m-%d %H:%M:%S] SamToFastq completed successfully")
    }

    output {
        File fastq1="${prefix}_1.fastq.gz"
        File fastq2="${prefix}_2.fastq.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "${memory}GB"
        disks: "local-disk ${diskGB} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}

# fastqc
task fastqc {

    File fastq1
    File? fastq2

    Float memory
    # Int disk_space
    Int num_threads
    Int num_preempt
    Int diskGB = ceil((size(fastq1) + size(fastq2)) * 1.25)

    String fastq1_name = sub(sub(basename(fastq1), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )
    String fastq2_name = sub(sub(basename(fastq2), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )

    command <<<
        set -euo pipefail
        echo $(date +"[%Y-%m-%d %H:%M:%S] Running FastQC...")
        fastqc ${fastq1} ${fastq2} \
            --threads ${num_threads} \
            --outdir .
        unzip -p ${fastq1_name}_fastqc.zip ${fastq1_name}_fastqc/fastqc_data.txt | gzip > ${fastq1_name}.fastqc_data.txt.gz
        unzip -p ${fastq2_name}_fastqc.zip ${fastq2_name}_fastqc/fastqc_data.txt | gzip > ${fastq2_name}.fastqc_data.txt.gz
        echo $(date +"[%Y-%m-%d %H:%M:%S] FastQC completed successfully")
    >>>

    output {
        File fastq1_fastqc_html = "${fastq1_name}_fastqc.html"
        File fastq1_fastqc_zip =  "${fastq1_name}_fastqc.zip"
        File fastq1_fastqc_data = "${fastq1_name}.fastqc_data.txt.gz"
        File fastq2_fastqc_html = "${fastq2_name}_fastqc.html"
        File fastq2_fastqc_zip =  "${fastq2_name}_fastqc.zip"
        File fastq2_fastqc_data = "${fastq2_name}.fastqc_data.txt.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "${memory}GB"
        disks: "local-disk ${diskGB} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}

# star
task star {
    File fastq1
    File? fastq2
    String prefix
    File star_index

    # STAR options
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

    Int memory
    # Int disk_space
    Int num_threads
    Int num_preempt
    Int diskGB = ceil((size(fastq1) + size(fastq2) + size(star_index)) * 2.5)

    command {
        set -euo pipefail

        if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]]; then
            tar -xvvf ${fastq1}
            fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            fastq2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            if [[ $fastq1_abs == *"*_1.fastq*" ]]; then  # no paired-end FASTQs found; check for single-end FASTQ
                fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
                fastq2_abs=''
            fi
        else
            # make sure paths are absolute
            fastq1_abs=${fastq1}
            fastq2_abs=${fastq2}
            if [[ $fastq1_abs != /* ]]; then
                fastq1_abs=$PWD/$fastq1_abs
                fastq2_abs=$PWD/$fastq2_abs
            fi
        fi

        echo "FASTQs:"
        echo $fastq1_abs
        echo $fastq2_abs

        # extract index
        echo $(date +"[%Y-%m-%d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1

        mkdir star_out
        # placeholders for optional outputs
        touch star_out/${prefix}.Aligned.toTranscriptome.out.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam.bai
        touch star_out/${prefix}.ReadsPerGene.out.tab  # run_STAR.py will gzip

        echo $(date +"[%Y-%m-%d %H:%M:%S] Running STAR alignment...")
        /src/run_STAR.py \
            star_index $fastq1_abs $fastq2_abs ${prefix} \
            --output_dir star_out \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \
            ${"--outFilterType " + outFilterType} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ${"--outSAMstrandField " + outSAMstrandField} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--quantMode " + quantMode} \
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--varVCFfile " + varVCFfile} \
            ${"--waspOutputMode " + waspOutputMode} \
            ${"--chimSegmentMin " + chimSegmentMin} \
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--chimOutJunctionFormat " + chimOutJunctionFormat} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            --threads ${num_threads}
        echo $(date +"[%Y-%m-%d %H:%M:%S] STAR alignment completed successfully")
    }

    output {
        File bam_file = "star_out/${prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "star_out/${prefix}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "star_out/${prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "star_out/${prefix}.Chimeric.out.junction.gz"
        File chimeric_bam_file = "star_out/${prefix}.Chimeric.out.sorted.bam"
        File chimeric_bam_index = "star_out/${prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "star_out/${prefix}.ReadsPerGene.out.tab.gz"
        File junctions = "star_out/${prefix}.SJ.out.tab.gz"
        File junctions_pass1 = "star_out/${prefix}._STARpass1/${prefix}.SJ.pass1.out.tab.gz"
        Array[File] logs = ["star_out/${prefix}.Log.final.out", "star_out/${prefix}.Log.out", "star_out/${prefix}.Log.progress.out"]
    }

    runtime {
        # V11 used 2.7.11a that doesn't have quantTranscriptomeSAMoutput
        # V11 uses 2.7.11b that does have quantTranscriptomeSAMoutput BanSingleEnd_ExtendSoftclip
        # https://github.com/alexdobin/STAR/releases
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11" 
        memory: "${memory}GB"
        disks: "local-disk ${diskGB} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}

# markduplicates
task markduplicates {

    File input_bam
    String prefix
    Int? max_records_in_ram
    Float? sorting_collection_size_ratio

    Float memory
    Int java_memory = floor(memory - 0.5)
    # Int disk_space
    Int num_threads
    Int num_preempt
    Int diskGB = ceil(size(input_bam) * 2.5)

    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")

    command {
        set -euo pipefail

        echo $(date +"[%Y-%m-%d %H:%M:%S] Running MarkDuplicates...")
        python3 -u /src/run_MarkDuplicates.py ${input_bam} ${prefix} \
            --memory ${java_memory} \
            ${"--max_records_in_ram " + max_records_in_ram} \
            ${"--sorting_collection_size_ratio " + sorting_collection_size_ratio}
        echo $(date +"[%Y-%m-%d %H:%M:%S] MarkDuplicates completed successfully")

        echo $(date +"[%Y-%m-%d %H:%M:%S] Indexing marked duplicate bam...")
        samtools index ${output_bam}
        echo $(date +"[%Y-%m-%d %H:%M:%S] Indexing completed successfully")
    }

    output {
        File bam_file = "${output_bam}"
        File bam_index = "${output_bam}.bai"
        File metrics = "${prefix}.marked_dup_metrics.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "${memory}GB"
        disks: "local-disk ${diskGB} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


# rsem
task rsem {

    File transcriptome_bam
    File rsem_reference
    String prefix

    Int memory
    # Int disk_space
    Int num_threads
    Int num_preempt
    Int diskGB = ceil(size(transcriptome_bam) * 1.25)

    Int? max_frag_len
    String? estimate_rspd
    String? is_stranded
    String? paired_end

    command {
        set -euo pipefail
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        echo $(date +"[%Y-%m-%d %H:%M:%S] Running RSEM calculate-expression...")
        /src/run_RSEM.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            ${"--paired_end " + paired_end} \
            --threads ${num_threads} \
            rsem_reference ${transcriptome_bam} ${prefix}
        gzip *.results
        echo $(date +"[%Y-%m-%d %H:%M:%S] RSEM calculate-expression completed successfully")
    }

    output {
        File genes="${prefix}.rsem.genes.results.gz"
        File isoforms="${prefix}.rsem.isoforms.results.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "${memory}GB"
        disks: "local-disk ${diskGB} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}

# rnaseqc2
task rnaseqc2 {

    File bam_file
    File genes_gtf
    String sample_id
    String? strandedness
    File? intervals_bed
    File? reference_fasta
    File? reference_fasta_index
    String? flags

    Int memory
    # Int disk_space
    Int num_threads
    Int num_preempt
    Int diskGB = ceil((size(bam_file) + size(genes_gtf) + size(intervals_bed) + size(reference_fasta)) * 1.5)

    command {
        set -euo pipefail
        echo $(date +"[%Y-%m-%d %H:%M:%S] Running RNA-SeQC 2...")
        touch ${sample_id}.fragmentSizes.txt
        touch ${sample_id}.gc_content.tsv
        rnaseqc ${genes_gtf} ${bam_file} . -s ${sample_id} ${"--bed " + intervals_bed} ${"--stranded " + strandedness} ${"--fasta " + reference_fasta} -vv ${flags}
        echo "  * compressing outputs"
        gzip *.gct
        echo $(date +"[%Y-%m-%d %H:%M:%S] RNA-SeQC 2 completed successfully")
    }

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
        File gene_counts = "${sample_id}.gene_reads.gct.gz"
        File exon_counts = "${sample_id}.exon_reads.gct.gz"
        File metrics = "${sample_id}.metrics.tsv"
        File gc_content = "${sample_id}.gc_content.tsv"
        File insertsize_distr = "${sample_id}.fragmentSizes.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11"
        memory: "${memory}GB"
        disks: "local-disk ${diskGB} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}