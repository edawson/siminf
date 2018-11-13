FROM ubuntu:16.04
MAINTAINER Eric T Dawson

RUN echo "deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse" | tee -a /etc/apt/sources.list && \
RUN apt-get update && apt-get install -yy gcc-4.9 wget tar bzip2 git python-dev build-essential zlib1g-dev vowpal-wabbit \
    tar git python-dev build-essential zlib1g-dev bzip2 gcc-4.9 libcurses

RUN mkdir /app
WORKDIR /app

RUN wget https://github.com/lomereiter/sambamba/releases/download/v0.6.5/sambamba_v0.6.5_linux.tar.bz2
RUN wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz

RUN tar xjf /app/sambamba_v0.6.5_linux.tar.bz2
RUN mv /app/sambamba_v0.6.5 /bin/sambamba

RUN tar xzf artbinmountrainier20160605linux64tgz.tgz
RUN mv /app/art_bin_MountRainier/art_* /bin/

RUN git clone --recursive https://github.com/edawson/siminf /app/siminf
RUN cp /app/siminf/scripts/* /bin/

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
RUN tar xzf bedtools-2.25.0.tar.gz && cd bedtools2 && make && cp bin/* /bin/

RUN git clone --recursive https://github.com/nh13/DWGSIM.git && cd DWGSIM && \
    make && cp dwgsim /bin/

RUN git clone --recursive https://github.com/edawson/rkmh.git && cd rkmh && make -j 3 && cp rkmh /bin/
ENV PATH="k/bin:${PATH}"
