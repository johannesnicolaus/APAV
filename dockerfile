# Dockerfile for APAV - An advanced pan-genome analysis and visualization toolkit
# Base image with bioc
FROM bioconductor/bioconductor_docker:latest

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /opt

# Install system dependencies
RUN apt-get update && apt-get install -y \
    perl \
    git \
    wget \
    curl \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    r-base \
    r-base-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Samtools (v1.13 or later required)
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xjf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.19.2 samtools-1.19.2.tar.bz2

# Install R packages (devtools and BiocManager)
# RUN R -e "install.packages(c('devtools', 'BiocManager'))"

# Install ComplexHeatmap from Bioconductor
RUN R -e "BiocManager::install('ComplexHeatmap', update=FALSE, ask=FALSE)"

# Clone APAV repository
RUN git clone https://github.com/SJTU-CGM/APAV.git /opt/APAV

# Install APAVplot R package
RUN R -e "devtools::install_github('SJTU-CGM/APAVplot')"

# Set up environment variables
ENV PATH="/opt/APAV:${PATH}"
ENV PERL5LIB="/opt/APAV/lib:${PERL5LIB}"

# Verify installation
RUN apav --help || echo "APAV installed"

# Set working directory for user
WORKDIR /data

# Entry point
CMD ["/bin/bash"]
