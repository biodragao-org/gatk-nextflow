#! /usr/bin/env nextflow

fastq_1 = file(params.fastqr1)
fastq_2 = file(params.fastqr2)

sample_name = params.sample_name
library_name = params.library_name
compression_level = params.compression_level
gatk_path = params.gatk_path
gotc_path = params.gotc_path
bwa_commandline = params.bwa_commandline
bwa_path = params.bwa_path

java_opt_convertubams = params.java_opt_convertubams
java_opt_bwamem = params.java_opt_bwamem
java_opt_mergebams = params.java_opt_mergebams
java_opt_markdups = params.java_opt_markdups
java_opt_sort = params.java_opt_sort
java_opt_fix = params.java_opt_fix
java_opt_baserecal = params.java_opt_baserecal
java_opt_bqsrreport = params.java_opt_bqsrreport
java_opt_applybqsr = params.java_opt_applybqsr
java_opt_gatherbams = params.java_opt_gatherbams
java_opt_haplotype = params.java_opt_haplotype
java_opt_mergegvcfs = params.java_opt_mergegvcfs

ref_fasta = file(params.fasta_ref)
ref_fasta_index = file( params.fasta_ref+'.fai' )
ref_alt = file( params.fasta_ref+'.64.alt' )
ref_amb = file( params.fasta_ref+'.64.amb' )
ref_ann = file( params.fasta_ref+'.64.ann' )
ref_bwt = file( params.fasta_ref+'.64.bwt' )
ref_pac = file( params.fasta_ref+'.64.pac' )
ref_sa = file( params.fasta_ref+'.64.sa' )
ref_dict = file( params.fasta_ref.replace(".fasta",".dict") )

dbSNP_vcf = file(params.dbsnp_vcf)
dbSNP_vcf_index = file(params.dbsnp_vcf+'.idx')
indels_mills_vcf = file(params.indels_mills)
indels_mills_vcf_index = file(params.indels_mills+'.tbi')
indels_humans_vcf = file(params.indels_humans)
indels_humans_vcf_index = file(params.indels_humans+'.tbi')

scattered_calling_intervals_list = file(params.scattered_calling_intervals_list)
sequence_grouping = file(params.sequence_grouping_interval_list)
sequence_grouping_with_unmapped = file(params.sequence_grouping_with_unmapped_interval_list)

outpath = params.outpath

process GetBwaVersion {

    output:
    stdout channel into outputs_bwa_version

    container "job-definition://gatk-bwaversion"

    """
    ${bwa_path}bwa 2>&1 | grep -e '^Version' | sed 's/Version: //' | tr -d '\n'
    """
}

process SplitFlowcellAndLaneR1 {

    input:
    file fastq_1

    output:
    file "*_R1.fq" into outputs_lane_r1

    container "job-definition://gatk-splitflowcellandlane"

    """
    zcat ${fastq_1} | awk  -v sample_name="$sample_name" 'BEGIN {FS = ":"} {flowcell=\$3; lane=\$4; print > sample_name"_"flowcell"."lane"_R1.fq"; for (i = 1; i <= 3; i++) {getline; print > sample_name"_"flowcell"."lane"_R1.fq"}}'
    """
}

process SplitFlowcellAndLaneR2 {

    input:
    file fastq_2

    output:
    file "*_R2.fq" into outputs_lane_r2

    container "job-definition://gatk-splitflowcellandlane"

    """
    zcat ${fastq_2} | awk  -v sample_name="$sample_name" 'BEGIN {FS = ":"} {flowcell=\$3; lane=\$4; print > sample_name"_"flowcell"."lane"_R2.fq"; for (i = 1; i <= 3; i++) {getline; print > sample_name"_"flowcell"."lane"_R2.fq"}}'
    """
}

process ConvertToUbams {

    input:
    file read1 from outputs_lane_r1.flatten()
    file read2 from outputs_lane_r2.flatten()

    output:
    stdout into readgroup_name
    file "${sample_name}_${read1.baseName.replaceAll(/${sample_name}_/, "").replaceAll(/_R1/, "")}.unmapped.bam" into outputs_unmapped_bam

    container "job-definition://gatk-convertunmappedbams"

    script:
    """
    RGN="${read1.baseName.replaceAll(/${sample_name}_/, "").replaceAll(/_R1/, "")}"
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_convertubams}" \
      FastqToSam \
      --FASTQ ${read1} \
      --FASTQ2 ${read2} \
      --OUTPUT ${sample_name}_${read1.baseName.replaceAll(/${sample_name}_/, "").replaceAll(/_R1/, "")}.unmapped.bam \
      --READ_GROUP_NAME ${read1.baseName.replaceAll(/${sample_name}_/, "").replaceAll(/_R1/, "")} \
      --SAMPLE_NAME ${sample_name} \
      --LIBRARY_NAME ${library_name} \
      --PLATFORM Illumina
    printf ${read1.baseName.replaceAll(/${sample_name}_/, "").replaceAll(/_R1/, "")}
    """
}

process UBamIndex {

    input:
    file umapped_bam from outputs_unmapped_bam

    output:
    file "${umapped_bam}.bai" into umapped_bam_index
    file "${umapped_bam}" into outputs_unmapped_bam_indexed

    container "job-definition://gatk-unmappedbamsindex"

    """
    samtools index -b ${umapped_bam} ${umapped_bam}.bai
    """
}

process SamToFastqAndBwaMem {

    input:
    file input_bam from outputs_unmapped_bam_indexed
    file input_bam_index from umapped_bam_index
    file ref_fasta
    file ref_fasta_index
    file ref_dict
    file ref_alt
    file ref_amb
    file ref_ann
    file ref_bwt
    file ref_pac
    file ref_sa

    output:
    file "${sample_name}.unmerged.bam" into output_unmerged_bam
    file "${sample_name}.bwa.stderr.log" into bwa_stderr_log
    file input_bam into output_unmapped_single_bam

    container "job-definition://gatk-bwa"

    script:
    """
    set -o pipefail
    set -e

    ls -la
	java -Dsamjdk.compression_level=${compression_level} ${java_opt_bwamem} -jar ${gotc_path}picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=/dev/stdout \
			INTERLEAVE=true \
			NON_PF=true \
    | \
		${bwa_path}bwa mem -K 100000000 -p -v 3 -t 16 -Y ${ref_fasta} /dev/stdin -  2> >(tee ${sample_name}.bwa.stderr.log >&2) \
    | \
		samtools view -1 - > ${sample_name}.unmerged.bam
    """
}

bamidgen = Channel.from( 1 .. 20 )
process MergeBamAlignment {

    input:
    file unmapped_bam from output_unmapped_single_bam
    val bwa_version from outputs_bwa_version
    file aligned_bam from output_unmerged_bam
    file ref_fasta
    file ref_fasta_index
    file ref_dict
    val fileid from bamidgen

    output:
    file "${sample_name}_${fileid}.aligned.unsorted.bam" into output_aligned_unsorted_bam

    container "job-definition://gatk-mergebam"

    script:
    """
    bash_fasta_ref=${ref_fasta}

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_mergebams}" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${aligned_bam} \
      --UNMAPPED_BAM ${unmapped_bam} \
      --OUTPUT ${sample_name}_${fileid}.aligned.unsorted.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER "unsorted" \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 2000000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t 16 -Y ${ref_fasta}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true
    """
}

process MarkDuplicates {

    input:
    file input_bams from output_aligned_unsorted_bam.collect(sort:true)

    output:
    file "${sample_name}_merged.bam" into output_bam_markeddups
    file "${sample_name}_merged_metrics.txt" into duplicate_metrics

    container "job-definition://gatk-markduplicates"

    script:
    """
    VARLIST=\$(echo "${input_bams}"|tr " " "\\n"|sort|tr "\\n" " "| sed -e 's/[[:space:]]*\$//')
    FILES=\${VARLIST//" "/" --INPUT "}

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_markdups}" \
      MarkDuplicates \
      --INPUT \$FILES \
      --OUTPUT ${sample_name}_merged.bam \
      --METRICS_FILE ${sample_name}_merged_metrics.txt \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
    """
}

process Sort {

    input:
    file input_bam from output_bam_markeddups
    file ref_dict
    file ref_fasta
    file ref_fasta_index

    output:
    file "${sample_name}_sorted.bam" into output_bam_sorted
    file "${sample_name}_sorted.bai" into output_bam_sorted_index
    file "${sample_name}_sorted.bam.md5" into output_bam_sorted_md5

    container "job-definition://gatk-sortandfixtags"

    script:
    """
    set -o pipefail

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_sort}" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT ${sample_name}_sorted.bam \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
    """
}

process FixTags {

    input:
    file input_bam from output_bam_sorted
    file input_bam_index from output_bam_sorted_index
    file input_bam_md5 from output_bam_sorted_md5
    file ref_dict
    file ref_fasta
    file ref_fasta_index

    output:
    file "${sample_name}_sorted_tagged.bam" into output_bam_sorted_tagged
    file "${sample_name}_sorted_tagged.bai" into output_bam_sorted_tagged_index
    file "${sample_name}_sorted_tagged.bam.md5" into output_bam_sorted_tagged_md5

    container "job-definition://gatk-sortandfixtags"

    """
    set -o pipefail

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_fix}" \
      SetNmAndUqTags \
      --INPUT ${input_bam} \
      --OUTPUT ${sample_name}_sorted_tagged.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
    """
}

idgen = Channel.from( 1 .. 20 )
process BaseRecalibrator {

    input:
    file input_bam from output_bam_sorted_tagged
    file input_bam_index from output_bam_sorted_tagged_index
    file dbSNP_vcf
    file dbSNP_vcf_index
    file indels_mills_vcf
    file indels_mills_vcf_index
    file indels_humans_vcf
    file indels_humans_vcf_index
    file ref_dict
    file ref_fasta
    file ref_fasta_index
    file sequence_group_interval from Channel.fromPath(sequence_grouping.splitText())
    val fileid from idgen

    output:
    file "${sample_name}_recalibration_report_${fileid.toString().padLeft(2, '0')}" into recalibration_report

    container "job-definition://gatk-baserecalibrator"

    """
    ${gatk_path} --java-options "${java_opt_baserecal}" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O "${sample_name}_recalibration_report_${fileid.toString().padLeft(2, '0')}" \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${indels_mills_vcf} \
      --known-sites ${indels_humans_vcf} \
      -L ${sequence_group_interval}
    """
}

process GatherBqsrReports {

    input:
    file input_bqsr_reports from recalibration_report.collect(sort:true)

    output:
    file "${sample_name}_bqsr_report" into output_bqsr_report

    container "job-definition://gatk-gatherbqsr"
    publishDir "${outpath}/BaseRecalibrator/", mode: 'copy', overwrite: true

    """
    ${gatk_path} --java-options "${java_opt_bqsrreport}" \
      GatherBQSRReports \
      -I ${input_bqsr_reports.join(" -I ")} \
      -O ${sample_name}_bqsr_report
    """
}

idgen2 = Channel.from( 1 .. 50 )
process ApplyBQSR {

    input:
    file input_bam from output_bam_sorted_tagged
    file input_bam_index from output_bam_sorted_tagged_index
    file recalibration_report from output_bqsr_report
    file sequence_group_interval from Channel.fromPath(sequence_grouping_with_unmapped.splitText())
    file ref_dict
    file ref_fasta
    file ref_fasta_index
    val fileid from idgen2

    output:
    file "${sample_name}_${fileid.toString().padLeft(2, '0')}.bam" into recalibrated_bam

    container "job-definition://gatk-applybqsr"

    """
    ${gatk_path} --java-options "${java_opt_applybqsr}" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${sample_name}_${fileid.toString().padLeft(2, '0')}.bam \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities \
      -L ${sequence_group_interval}
    """
}

process GatherBamFiles {

    input:
    file input_bams from recalibrated_bam.collect(sort:true)

    output:
    file "${sample_name}.bam" into output_bam
    file "${sample_name}.bai" into output_bam_index
    file "${sample_name}.bam.md5" into output_bam_md5

    container "job-definition://gatk-gatherbams"
    publishDir "${outpath}/bam/", mode: 'copy', overwrite: true

    """
    VARLIST=\$(echo "${input_bams}"|tr " " "\\n"|sort|tr "\\n" " "| sed -e 's/[[:space:]]*\$//')
    FILES=\${VARLIST//" "/" --INPUT "}
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_gatherbams}" \
      GatherBamFiles \
      --INPUT \$FILES \
      --OUTPUT ${sample_name}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
    """
}

idgen3 = Channel.from( 1 .. 50 )
process HaplotypeCaller {

    input:
    file input_bam from output_bam
    file input_bam_index from output_bam_index
    file ref_fasta
    file ref_fasta_index
    file ref_dict
    file interval_list from Channel.fromPath(scattered_calling_intervals_list.splitText())
    val fileid from idgen3

    output:
    file "${sample_name}_${fileid.toString().padLeft(2, '0')}.g.vcf" into output_vcf
    file "${sample_name}_${fileid.toString().padLeft(2, '0')}.g.vcf.idx" into output_vcf_index

    container "job-definition://gatk-haplotype"

    """
    set -e

    ${gatk_path} --java-options "${java_opt_haplotype}" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --output ${sample_name}_${fileid.toString().padLeft(2, '0')}.g.vcf \
      -contamination 0 \
      -ERC GVCF \
      -L ${interval_list}
    """
}

process MergeGVCFs {

    input:
    file input_vcfs from output_vcf.collect(sort:true)
    file input_vcfs_indexes from output_vcf_index.collect(sort:true)

    output:
    file "${sample_name}.g.vcf" into output_final_vcf
    file "${sample_name}.g.vcf.idx" into output_final_vcf_index

    container "job-definition://gatk-mergegvcfs"
    publishDir "${outpath}/gVCF/", mode: 'copy', overwrite: true

    """
    set -e

    VARLIST=\$(echo "${input_vcfs}"|tr " " "\\n"|sort|tr "\\n" " "| sed -e 's/[[:space:]]*\$//')
    FILES=\${VARLIST//" "/" --INPUT "}
    ${gatk_path} --java-options "${java_opt_mergegvcfs}"  \
      MergeVcfs \
      --INPUT \$FILES \
      --OUTPUT ${sample_name}.g.vcf
    """
}
