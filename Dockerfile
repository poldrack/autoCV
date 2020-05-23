# Dockerfile for ThinkStats

FROM python:3.6-stretch

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
    pypatent

WORKDIR /data
CMD ["/usr/bin/xelatex", "autocv_template.tex"]
