FROM ubuntu:16.04
MAINTAINER Eric T Dawson

RUN apt-get update && apt-get install -yy gcc-4.9 wget tar bzip2 git python-dev build-essential zlib1g-dev

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

