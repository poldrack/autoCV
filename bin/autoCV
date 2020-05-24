#!/bin/bash

# generate CV files
docker run -it -v $PWD:/data  poldrack/latex-python python /usr/local/bin/make_cv.py

# render C
docker run -it -v $PWD:/data  poldrack/latex-python

