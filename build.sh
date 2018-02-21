#!/bin/bash

NAME=thesis_main
FINALNAME=thomasstainerthesis2015

pdflatex ${NAME}.tex
pdflatex ${NAME}.tex
bibtex ${NAME}.aux
pdflatex ${NAME}.tex
mv ${NAME}.pdf ${FINALNAME}.pdf

# cleanup
rm ${NAME}.aux ${NAME}.bbl ${NAME}.blg ${NAME}.idx ${NAME}.log ${NAME}.nlo ${NAME}.out ${NAME}.toc
