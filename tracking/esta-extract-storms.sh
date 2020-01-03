#!/bin/bash

# Script to A) track storms and then B) extract snowfall/calculate RSI values and other statistics
# from LENS member (NUMID). 
#
# Colin Zarzycki (zarzycki@ucar.edu)
#
# Usage is ./lens-snow-tracker $NUMID
# Driven by "batch-loginnodes.sh" for parallel processing of LENS data

starttime=$(date -u +"%s")

################################################################################

NAMELIST=$1

while IFS='=' read -r var value ; do
  value=`echo $value | sed 's/\(.*\),/\1/'`  #strip off last comma
  value=`echo $value | sed 's/\"//g'`
  echo "... SETTING: ${var} to ${value}"
  export "$var"="$value"
done < ${NAMELIST}

export ESTAPATH

################################################################################

ncl ${ESTAPATH}/tracking/extract_individual_storms-lite.ncl 'nlfile="'${NAMELIST}'"'

################################################################################

# Calculate script timing
endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "`basename $0` ${TE_UQSTR} ${tottime}\n" >> timing.txt
