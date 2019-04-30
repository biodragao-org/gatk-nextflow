#!/usr/bin/env bash
pwd

STACKNAME="gatk-DNA-seq"
TEMPLATE="https://s3.amazonaws.com/gfb-genomics/gatk_dnaseq/cfn/iam.template.yaml"

aws cloudformation create-stack --stack-name "${STACKNAME}" --template-url "${TEMPLATE}" --capabilities CAPABILITY_NAMED_IAM
