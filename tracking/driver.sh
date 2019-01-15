#!/bin/bash

module load parallel

TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

declare -a ENSMEMS=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20")

YYYY=2019
MM=01
DD=14
HH=12
YYYYMMDDHH=${YYYY}${MM}${DD}${HH}

rm ${COMMANDFILE}
for zz in "${ENSMEMS[@]}"; do
  THECOMMAND="/glade/u/home/zarzycki/proc-GEFS/batch.sh ${zz} ${YYYY} ${MM} ${DD} ${HH}"
  echo ${THECOMMAND} >> ${COMMANDFILE}
done

parallel --jobs 3 -u < ${COMMANDFILE}

rm ${COMMANDFILE}
for zz in "${ENSMEMS[@]}"; do
  cp nl.gefs.base nl.gefs.${zz}
  sed -i "s|_ENSMEMBER_|${zz}|g" nl.gefs.${zz}
  NAMELIST=nl.gefs.${zz}
  THECOMMAND="./esta-gen-filelist.sh ${NAMELIST} ; ./esta-etc-tracking.sh ${NAMELIST} ; ./esta-extract-storms.sh ${NAMELIST} ; ./esta-calc-RSI.sh ${NAMELIST}"
  echo ${THECOMMAND} >> ${COMMANDFILE}
done

parallel --jobs 3 -u < ${COMMANDFILE}

cd /glade/u/home/zarzycki/proc-GEFS/
rm highestRSI.txt
for zz in "${ENSMEMS[@]}"; do
  ncl select-highest-RSI.ncl 'ENSMEM="'${zz}'"'
done

ncl plot-highest-RSI-dist.ncl 'DATE="'${YYYYMMDDHH}'"'

### CLEANUP
rm -rf /glade/u/home/zarzycki/scratch/GEFS/*

cd /glade/u/home/zarzycki/snow-tracking/tracking/
rm traj.GEFS*
rm filelist.PSL.GEFS.*
rm nl.gefs.??
rm ${COMMANDFILE}