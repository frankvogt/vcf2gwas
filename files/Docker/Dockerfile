FROM python:slim
WORKDIR /vcf2gwas
RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && echo "Running $(conda --version)" && \
    conda init bash && \
    . /root/.bashrc && \
    conda update conda && \
    conda create -n python-app && \
    conda activate python-app && \
    conda install vcf2gwas -c bioconda -c conda-forge -c fvogt257 && \
    conda list && \
    vcf2gwas --version
COPY ./commands.sh /
RUN chmod +x /commands.sh
ENTRYPOINT ["/commands.sh"]
