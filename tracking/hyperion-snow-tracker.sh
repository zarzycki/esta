#!/bin/bash

# Script to A) track storms and then B) extract snowfall/calculate RSI values and other statistics
# from LENS member (NUMID). 
#
# Colin Zarzycki (zarzycki@ucar.edu)
#
# Usage is ./lens-snow-tracker $NUMID
# Driven by "batch-loginnodes.sh" for parallel processing of LENS data

date

export ESTAPATH="/global/homes/c/czarzyck/esta/"

SWE="12"
DESCSTR=HYPERION
RSIOUTDIR=/global/homes/c/czarzyck/scratch/hyperion/
TRAJFILE="./traj"
RSISNOWFILE="RSI.SNOW."${DESCSTR}".csv"

TRAJDIR=${ESTAPATH}/tracking/
RSIDIR=${ESTAPATH}/calc_RSI/

starttime=$(date -u +"%s")

###################################################################################################

starttimeextract=$(date -u +"%s")
EXTRACTOUTFILE=${EXTRACTOUTFILE}".tempest.nc"
ncl ${TRAJDIR}/extract_individual_storms-lite.ncl
endtimeextract=$(date -u +"%s")

###################################################################################################

starttimeRSI=$(date -u +"%s")
mkdir -p ${RSIOUTDIR}
SWES=(${SWE} ${SWE})
SNOWVARNAMES=("SUM_PRECB_SN" "NONE")
RSISNOWFILES=(${RSISNOWFILE}".SNOW.csv" ${RSISNOWFILE}".PRECT.csv")

# LOOPVAR does both snow and total precip (pretending total precip is snow)
for LOOPVAR in `seq 0 1`
do
  rm ${RSIOUTDIR}/${RSISNOWFILES[${LOOPVAR}]}
  LOOPST=0
  LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
  for IX in `seq ${LOOPST} ${LOOPEN}`
  do
     echo $IX
    (set -x; ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} SWE=${SWES[${LOOPVAR}]} \
      'imgDir="'${RSIOUTDIR}'/images-RSI-tempest/'${DESCSTR}'/"' \
      'RSIoutFile="'${RSISNOWFILES[${LOOPVAR}]}'"' \
      'SNOWVARNAME="'${SNOWVARNAMES[${LOOPVAR}]}'"')
  done
done
endtimeRSI=$(date -u +"%s")