import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:fastqc_v1-0_BETA/versions/2/plain-WDL/descriptor" as fastqc_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:star_v1-0_BETA/versions/8/plain-WDL/descriptor" as star_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:markduplicates_v1-0_BETA/versions/6/plain-WDL/descriptor" as markduplicates_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rsem_v1-0_BETA/versions/6/plain-WDL/descriptor" as rsem_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rnaseqc2_v1-0_BETA/versions/4/plain-WDL/descriptor" as rnaseqc_wdl

#---------------------------------------------------------------------------------------------------------------
# TASKS 
#---------------------------------------------------------------------------------------------------------------

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
    Int disk_space
    Int num_threads
    Int num_preempt

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
        echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1

        mkdir star_out
        # placeholders for optional outputs
        touch star_out/${prefix}.Aligned.toTranscriptome.out.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam.bai
        touch star_out/${prefix}.ReadsPerGene.out.tab  # run_STAR.py will gzip

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
        # V10 used 2.7.11a that doesn't have quantTranscriptomeSAMoutput
        # V11 uses 2.7.11b that does have quantTranscriptomeSAMoutput BanSingleEnd_ExtendSoftclip
        # https://github.com/alexdobin/STAR/releases
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V11" 
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}

#---------------------------------------------------------------------------------------------------------------
# WORKFLOW 
#---------------------------------------------------------------------------------------------------------------
workflow rnaseq_pipeline_fastq_workflow {

    String prefix

    call fastqc_wdl.fastqc {}

    call star {
        input: prefix=prefix
    }

    call markduplicates_wdl.markduplicates {
        input: input_bam=star.bam_file, prefix=prefix
    }

    call rsem_wdl.rsem {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call rnaseqc_wdl.rnaseqc2 {
        input: bam_file=markduplicates.bam_file, sample_id=prefix
    }
}
