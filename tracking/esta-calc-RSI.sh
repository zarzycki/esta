#!/bin/bash

# Script to A) track storms and then B) extract snowfall/calculate RSI values and other statistics
# from LENS member (NUMID). 
#
# Colin Zarzycki (zarzycki@ucar.edu)
#
# Usage is ./lens-snow-tracker $NUMID
# Driven by "batch-loginnodes.sh" for parallel processing of LENS data

date

starttime=$(date -u +"%s")

################################################################################

NAMELIST=$1

while IFS='=' read -r var value ; do
  value=`echo $value | sed 's/\(.*\),/\1/'`  #strip off last comma
  value=`echo $value | sed 's/\"//g'`
  echo "... SETTING: ${var} to ${value}"
  export "$var"="$value"
done < ${NAMELIST}

################################################################################

# Hardcoded options
RSISNOWFILE="RSI."${TE_UQSTR}
RSIDIR=${ESTAPATH}/calc_RSI/
TRAJFILE=${TRAJDIR}"/traj."${TE_UQSTR}

################################################################################

mkdir -p ${RSI_OUTDIR}
SWES=(${SWE} ${SWE})
SNOWVARNAMES=($RSI_SUMVAR "NONE")
RSISNOWFILES=(${RSISNOWFILE}".SNOW.csv" ${RSISNOWFILE}".PRECT.csv")

# LOOPVAR does snow first, if you'd like to do another variable (like PRECT/RPI from
# Zarzycki (2018), GRL, then change above values and change to seq 0 1
for LOOPVAR in `seq 0 0`
do
  rm ${RSI_OUTDIR}/${RSISNOWFILES[${LOOPVAR}]}
  LOOPST=0
  LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
  for IX in `seq ${LOOPST} ${LOOPEN}`
  do
     echo $IX
    (set -x; ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} SWE=${SWES[${LOOPVAR}]} \
      'nlfile="'${NAMELIST}'"' \
      'imgDir="'${RSI_OUTDIR}'/images-RSI-tempest/'${TE_UQSTR}'/"' \
      'RSIoutFile="'${RSISNOWFILES[${LOOPVAR}]}'"' \
      'SNOWVARNAME="'${SNOWVARNAMES[${LOOPVAR}]}'"')
  done
done

################################################################################

# Calculate script timing
endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "`basename $0` ${TE_UQSTR} ${tottime}\n" >> timing.txt
