# Dockerfile for autoCV

FROM python:3.8-buster

# apt-get installs
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    vim \
    wget \
    pandoc \
    texlive-full \
    python-pygments gnuplot \
    make git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# pip installs

RUN pip install \
    pandas \
    biopython \
    requests \
    crossrefapi \
    scholarly \
    pypatent \
    pytest \
    pytest-cov \
    flake8

## this forces rebuild each time, when build arg is set to date
ARG DUMMY=unknown
RUN DUMMY=${DUMMY} pip install git+https://github.com/poldrack/autoCV

WORKDIR /data
CMD ["/usr/bin/xelatex", "autocv_template.tex"]
