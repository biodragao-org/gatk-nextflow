# gatk-nextflow
A Nextflow workflow to run GATK 4.1 modeled after the Broad's WDLs

# Setup

## Create s3 buckets

```bash
aws s3 mb s3://my-workfdir
aws s3 mb s3://my-refdata
```

## Build docker images

### Nextflow

#### Login to ECR
```bash
$(aws ecr get-login --no-include-email --region us-east-1)
```

#### Build the image and push it to the repository

Update the ECR repository reference below

```bash
cd docker/nextflow
export REPO=123456789.dkr.ecr.us-east-1.amazonaws.com/nextflow
export VERSION=19.04.0

docker build -t ${REPO}:${VERSION} --build-arg VERSION=${VERSION} --rm .
docker push ${REPO}:${VERSION}
cd ../..
```

### Genomes in the Cloud

#### Login to ECR
```bash
$(aws ecr get-login --no-include-email --region us-east-1)
```

#### Build the image and push it to the repository

Update the ECR repository reference below

```bash
cd docker/genomes_in_the_cloud
export REPO=123456789.dkr.ecr.us-east-1.amazonaws.com/genomes-in-the-cloud
export VERSION=1.0.0

docker build -t ${REPO}:${VERSION} --rm .
docker push ${REPO}:${VERSION}
cd ../..
```

### Nextflow

## Copy reference data (hg38)

```bash
mkdir hg38
cd hg38/
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0 .
aws s3 sync . s3://my-refdata/broad-references/hg38/
cd ..
rm -fr hg38
```

## Build the CloudFormation stack

Update the ECR repository reference in iam.template.yaml prior to doing the following:

```bash
cd gatk_dnaseq
aws s3 cp iam.template.yaml s3://my-refdata/cfn/
./build_environment.sh
cd ..
```

# Generate Pre-requisites

## Run the workflow to produce some of the interval files required for the full workflow

```bash
cd gatk_dnaseq_pre
./run_pipeline.sh
cd ..
```

This workflow produces the following files:
```
s3://my-refdata/broad-references/hg38/sequence_groups/sequence_grouping_*.interval_list
s3://my-refdata/broad-references/hg38/sequence_groups/sequence_group_list.txt
s3://my-refdata/broad-references/hg38/sequence_groups_unmapped/sequence_grouping_with_unmapped_*.interval_list
s3://my-refdata/broad-references/hg38/sequence_groups_unmapped/sequence_grouping_with_unmapped_list.txt
```

## Sync the scattered interval list maintained by the Broad to s3

```bash
cd gatk_dnaseq
./get_interval_list.sh
aws s3 cp scattered_interval_list.txt s3://my-refdata/broad-references/hg38/interval_list/
```

# Run Workflow

There are two versions of the GATK workflow.  
The bwa_start does not convert the pair of FASTQ files into ubams and passes the pair of files directly into BWA.
The ubams does convert the pair of FASTQ files into ubams and more closely follows the Broad's best practices.
Post alignment, the two workflows are identical and are closely modeled after the Broad's WDLs best practices.

```bash
cd gatk_dnaseq/ubams
./run_pipeline.sh
```