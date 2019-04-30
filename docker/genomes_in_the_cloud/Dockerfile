FROM r-base
LABEL maintainer "Adam Tebbe <atebbe@goldfinchbio.com>"

ARG PICARD_VERSION=2.18.23
ARG BWA_VERSION=0.7.17
ARG SAMTOOLS_VERSION=1.9
ARG SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download"
ARG VERIFYBAMID_VERSION=1.0.5
ARG GATK_VERSION=4.1.0.0

RUN apt-get update && \
    apt-get remove -y openjdk-9-jre-headless openjdk-9-jdk-headless && \
    apt-get install -y python3 python3-pip openjdk-8-jre python3-boto3 build-essential openjdk-8-jre-headless openjdk-8-jdk wget libterm-readline-perl-perl libcurl4-openssl-dev libssl1.0.2 libssl-dev libgomp1 cmake make gcc g++ git libz-dev libbz2-dev automake autoconf awscli && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --upgrade awscli

RUN cd /root && \
    wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 && \
    bunzip2 bwa-${BWA_VERSION}.tar.bz2 && \
    tar -xvf bwa-${BWA_VERSION}.tar && \
    cd bwa-${BWA_VERSION} && \
    make && \
    mv bwa /usr/local/bin/ && \
    cd .. && \
    rm -fr bwa-${BWA_VERSION}

RUN cd /root && \
    wget --no-check-certificate ${SAMTOOLS_URL}/$SAMTOOLS_VERSION/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    make && \
    make install && \
    cd .. && \
    rm -fr samtools-${SAMTOOLS_VERSION}

RUN cd /root && \
    wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip && \
    mv gatk-${GATK_VERSION} /usr/local/bin/gatk && \
    rm gatk-${GATK_VERSION}.zip

RUN cd /root && \
    wget https://github.com/Griffan/VerifyBamID/archive/${VERIFYBAMID_VERSION}.tar.gz && \
    tar -xvzf ${VERIFYBAMID_VERSION}.tar.gz && \
    cd VerifyBamID-${VERIFYBAMID_VERSION} && \
    SHLIB_VERSION_NUMBER="1.0.2" cmake . && \
    make && \
    ls -l && \
    ls -l bin/ && \
    mv bin/* /usr/local/bin/ && \
    chmod 755 /usr/local/bin/VerifyBamID && \
    chmod 755 /usr/local/bin/run.plot.sh && \
    cd /root && \
    rm -fr VerifyBamID-${VERIFYBAMID_VERSION}

RUN pip3 install boto3 psutil
RUN wget -O /usr/local/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar
