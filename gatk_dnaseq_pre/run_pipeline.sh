#!/usr/bin/env bash
set -e

function runPipeline()
{
  BASEOUTPUT=s3://my-refdata/broad-references/hg38
  WORKDIR=s3://my-workdir/low_coverage_variants/gatk

  docker run --rm -v `pwd`:/mnt/data --workdir /mnt/data -v ~/.aws:/root/.aws 123456789.dkr.ecr.us-east-1.amazonaws.com/nextflow:19.04.0 \
  /usr/local/bin/nextflow run /mnt/data/prepare_workflow.nf \
  --outpath ${BASEOUTPUT} \
  -bucket-dir ${WORKDIR}
}

runPipeline
