#!/usr/bin/env bash
set -e

function runPipeline()
{
  SAMPLE=$1
  FASTQ1=$2
  FASTQ2=$3

  BASEOUTPUT=s3://my-outdir/
  WORKDIR=s3://my-workdir/

  docker run --rm -v `pwd`:/mnt/data --workdir /mnt/data -v ~/.aws:/root/.aws -e TMPDIR=/tmp 227114915345.dkr.ecr.us-east-1.amazonaws.com/nextflow:19.04.0 \
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

runPipeline GRS-0026552 s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.1333_R1.fastq.gz s3://gfb-registry-raw/controls/nist/human/dna/fastq/genome_in_a_bottle/illumina300x_.1333_R2.fastq.gz &
