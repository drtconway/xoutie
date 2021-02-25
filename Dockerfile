FROM ubuntu:groovy

ENV DEBIAN_FRONTEND="noninteractive"

RUN apt update && \
    apt install -y tzdata && \
    apt install -y \
        python3-pip \
        python3-setuptools \
        pypy3 \
        pypy3-dev \
        gfortran \
        liblapack-dev \
        libopenblas-dev \
        samtools && \
    pypy3 -m pip install --upgrade pip && \
    pypy3 -m pip install --upgrade setuptools

ADD xoutie-requirements.txt /tmp/
RUN pypy3 -m pip install -r /tmp/xoutie-requirements.txt && rm /tmp/xoutie-requirements.txt

ADD . /source/

RUN cd /source/ && pypy3 -m pip install . && cd / && rm -rf /source/
