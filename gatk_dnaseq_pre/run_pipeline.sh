#!/usr/bin/env bash
set -e

function runPipeline()
{
  SAMPLE=$1
  FASTQ1=$2
  FASTQ2=$3
  BASEOUTPUT=s3://gfb-genomics/broad-references/hg38
  WORKDIR=s3://gfb-registry-workdir/low_coverage_variants/gatk

  time docker run --rm -v `pwd`:/mnt/data --workdir /mnt/data -v ~/.aws:/root/.aws 227114915345.dkr.ecr.us-east-1.amazonaws.com/nextflow:19.01.0 \
  /usr/local/bin/nextflow run /mnt/data/prepare_workflow.nf \
  -with-report ${SAMPLE}_report.html \
  -with-trace ${SAMPLE}_trace.txt \
  -with-timeline ${SAMPLE}_timeline.html \
  -with-dag ${SAMPLE}_flowchart.png \
  --outpath ${BASEOUTPUT} \
  -bucket-dir ${WORKDIR}/${SAMPLE}/workdir \
  --output ${WORKDIR}/${SAMPLE}/output \
  --resume

  aws s3 cp ${SAMPLE}_report.html ${BASEOUTPUT}/${SAMPLE}/
  aws s3 cp ${SAMPLE}_trace.txt ${BASEOUTPUT}/${SAMPLE}/
  aws s3 cp ${SAMPLE}_timeline.html ${BASEOUTPUT}/${SAMPLE}/
  aws s3 cp ${SAMPLE}_flowchart.png ${BASEOUTPUT}/${SAMPLE}/
}

runPipeline
