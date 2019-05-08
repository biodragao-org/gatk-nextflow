#!/usr/bin/env bash

STACKNAME="gatk-DNA-seq"
TEMPLATE="https://s3.amazonaws.com/my-refdata/cfn/iam.template.yaml"

aws cloudformation create-stack --stack-name "${STACKNAME}" --template-url "${TEMPLATE}" --capabilities CAPABILITY_NAMED_IAM
