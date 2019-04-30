FROM openjdk:8-jre
LABEL maintainer "Adam Tebbe <atebbe@goldfinchbio.com>"

ARG VERSION="19.04.0"

RUN apt-get update && \
    apt-get install -y wget awscli graphviz apt-transport-https ca-certificates curl software-properties-common sendmail sendmail-cf m4 && \
    curl -fsSL https://download.docker.com/linux/debian/gpg | apt-key add - && \
    add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable" && \
    apt-get update && \
    apt-get install -y docker-ce && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN cd /root && \
    wget https://github.com/nextflow-io/nextflow/releases/download/v${VERSION}/nextflow-${VERSION}-all && \
    mv nextflow-${VERSION}-all /usr/local/bin/nextflow && \
    chmod 755 /usr/local/bin/nextflow

VOLUME /mnt/data
