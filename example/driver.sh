#!/bin/bash

ESTAPATH=/glade/home/zarzycki/snow-tracking/
NAMELIST=nl.jra.sample

${ESTAPATH}/tracking/esta-etc-tracking.sh ${NAMELIST} 
${ESTAPATH}/tracking/esta-extract-storms.sh ${NAMELIST}
${ESTAPATH}/tracking/esta-calc-RSI.sh ${NAMELIST}"