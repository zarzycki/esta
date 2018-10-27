#!/bin/bash

# Script to A) track storms and then B) extract snowfall/calculate RSI values and other statistics
# from LENS member (NUMID). 
#
# Colin Zarzycki (zarzycki@ucar.edu)
#
# Usage is ./lens-snow-tracker $NUMID
# Driven by "batch-loginnodes.sh" for parallel processing of LENS data

date

NAMELIST="nl.hyperion"

while IFS='=' read -r var value ; do
  value=`echo $value | sed 's/\(.*\),/\1/'`  #strip off last comma
  value=`echo $value | sed 's/\"//g'`
  echo "... SETTING: ${var} to ${value}"
  export "$var"="$value"
done < ${NAMELIST}

# Hardcoded options
RSISNOWFILE="RSI.SNOW."${DESCSTR}".csv"
TRAJDIR=${ESTAPATH}/tracking/
RSIDIR=${ESTAPATH}/calc_RSI/

starttime=$(date -u +"%s")

###################################################################################################

starttimeextract=$(date -u +"%s")
#ncl ${TRAJDIR}/extract_individual_storms-lite.ncl 'nlfile="'${NAMELIST}'"'
endtimeextract=$(date -u +"%s")

###################################################################################################

starttimeRSI=$(date -u +"%s")
mkdir -p ${RSI_OUTDIR}
SWES=(${SWE} ${SWE})
SNOWVARNAMES=("SUM_PRECB_SN" "NONE")
RSISNOWFILES=(${RSISNOWFILE}".SNOW.csv" ${RSISNOWFILE}".PRECT.csv")

# LOOPVAR does both snow and total precip (pretending total precip is snow)
for LOOPVAR in `seq 0 1`
do
  rm ${RSI_OUTDIR}/${RSISNOWFILES[${LOOPVAR}]}
  LOOPST=0
  LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
  for IX in `seq ${LOOPST} ${LOOPEN}`
  do
     echo $IX
    (set -x; ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} SWE=${SWES[${LOOPVAR}]} \
      'nlfile="'${NAMELIST}'"' \
      'imgDir="'${RSI_OUTDIR}'/images-RSI-tempest/'${DESCSTR}'/"' \
      'RSIoutFile="'${RSISNOWFILES[${LOOPVAR}]}'"' \
      'SNOWVARNAME="'${SNOWVARNAMES[${LOOPVAR}]}'"')
  done
done
endtimeRSI=$(date -u +"%s")
