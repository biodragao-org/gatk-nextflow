#!/usr/bin/env bash
set -e

function runPipeline()
{
  SAMPLE=$1
  FASTQ1=$2
  FASTQ2=$3
  BASEOUTPUT=s3://gfb-genomics/low_coverage_variants/gatk_try2
  WORKDIR=s3://gfb-registry-workdir/low_coverage_variants/gatk_try2

  time docker run --rm -v `pwd`:/mnt/data --workdir /mnt/data -v ~/.aws:/root/.aws -e TMPDIR=/tmp 227114915345.dkr.ecr.us-east-1.amazonaws.com/nextflow:19.01.0 \
  /usr/local/bin/nextflow run /mnt/data/processing-for-variant-discovery-gatk4.nf \
  -with-report ${SAMPLE}_report.html \
  -with-trace ${SAMPLE}_trace.txt \
  -with-timeline ${SAMPLE}_timeline.html \
  -with-dag ${SAMPLE}_flowchart.png \
  --outpath ${BASEOUTPUT}/ \
  -bucket-dir ${WORKDIR}/${SAMPLE}/workdir \
  --output ${WORKDIR}/${SAMPLE}/output \
  --fastqr1 ${FASTQ1} \
  --fastqr2 ${FASTQ2} \
  --sample_name ${SAMPLE}

  aws s3 cp ${SAMPLE}_report.html ${BASEOUTPUT}/${SAMPLE}/
  aws s3 cp ${SAMPLE}_trace.txt ${BASEOUTPUT}/${SAMPLE}/
  aws s3 cp ${SAMPLE}_timeline.html ${BASEOUTPUT}/${SAMPLE}/
  aws s3 cp ${SAMPLE}_flowchart.png ${BASEOUTPUT}/${SAMPLE}/
}

#sleep 10; runPipeline GRS-0000047 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000047_S2_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000047_S2_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000079 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000079_S1_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000079_S1_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000091 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000091_S4_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000091_S4_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000282 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000282_S8_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000282_S8_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000366 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000366_S4_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000366_S4_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000376 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000376_S5_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000376_S5_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000394 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000394_S6_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000394_S6_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000404 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000404_S7_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000404_S7_R2_001.fastq.gz &
sleep 10; runPipeline GRS-0000563 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000563_S5_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000563_S5_R2_001.fastq.gz &
sleep 10; runPipeline GRS-0000592 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000592_S7_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000592_S7_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000002 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000002_S3_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000002_S3_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000008 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000008_S4_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000008_S4_R2_001.fastq.gz &
sleep 10; runPipeline GRS-0000015 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000015_S6_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000015_S6_R2_001.fastq.gz &
sleep 10; runPipeline GRS-0000017 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000017_S4_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000017_S4_R2_001.fastq.gz &
sleep 10; runPipeline GRS-0000019 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000019_S8_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000019_S8_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000325 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000325_S1_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20171006_batch2/Sample_GRS0000325_S1_R2_001.fastq.gz &
#sleep 10; runPipeline GRS-0000013 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000013_S6_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000013_S6_R2_001.fastq.gz &
sleep 10; runPipeline GRS-0000136 s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000136_S5_R1_001.fastq.gz s3://gfb-registry-raw/samples/mount_sinai/human/dna/fastq/20170626_batch1/Sample_GRS0000136_S5_R2_001.fastq.gz &

sleep 10; runPipeline GRS-0026550 s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.00333_R1.fastq.gz s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.00333_R2.fastq.gz &
sleep 10; runPipeline GRS-0026551 s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.03333_R1.fastq.gz s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.03333_R2.fastq.gz &
#sleep 10; runPipeline GRS-0026552 s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.1333_R1.fastq.gz s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.1333_R2.fastq.gz &

wait
