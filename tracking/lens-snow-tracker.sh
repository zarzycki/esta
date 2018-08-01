#!/bin/bash

# Script to A) track storms and then B) extract snowfall/calculate RSI values and other statistics
# from LENS member (NUMID). 
#
# Colin Zarzycki (zarzycki@ucar.edu)
#
# Usage is ./lens-snow-tracker $NUMID
# Driven by "batch-loginnodes.sh" for parallel processing of LENS data

date

NUMID=$1

DO_TEMPEST=false           # Do we do Lagrangian/Tempest trajectories?
DO_CONNTRACK=false        # Do we do "connected" trajectories?
PROCESS_TEMPEST=false      # Do we extract/calc RSI for Lagrangian/Tempest trajectories?
PROCESS_CONNTRACK=false   # Do we extract/calc RSI for "connected" trajectories?
SEASONFRAC_TEMPEST=true
SEASONFRAC_CONNTRACK=true

SWE="12"                   # snow water equivalent (12 = default, 9-13 is probably acceptable range)
THRESH="5e-9"              # Connected trajectory threshold (area average precipitation)

MODEL="LENS"
#RUNDATES="1990010100Z-2005123118Z"
#RUNDATES="2026010100Z-2035123118Z"
RUNDATES="2071010100Z-2080123118Z"
#RUNCONFIG="B20TRC5CNBDRD"
RUNCONFIG="BRCP85C5CNBDRD"
CONFIG=${THRESH}"_"${SWE}

# Logic for LENS members done at U Toronto
if   [ ${NUMID} == '036' ]; then
  NUMID='101'
elif [ ${NUMID} == '037' ]; then
  NUMID='102'
elif [ ${NUMID} == '038' ]; then
  NUMID='103'
elif [ ${NUMID} == '039' ]; then
  NUMID='104'
elif [ ${NUMID} == '040' ]; then
  NUMID='105'
elif [ ${NUMID} == '041' ]; then
  NUMID='106'
elif [ ${NUMID} == '042' ]; then
  NUMID='107'
else
  NUMID=$NUMID
fi

DESCSTR=${MODEL}"."${RUNDATES:0:4}"."${NUMID}"."${CONFIG}

RSIDIR=/glade/u/home/zarzycki/snow-tracking/calc_RSI/
RSISNOWFILE="RSI.SNOW."${DESCSTR}".csv"
RSIPRECFILE="RSI.PREC."${DESCSTR}".csv"

NCLPATH=/glade/u/home/zarzycki/snow-tracking/tracking/
RSIOUTDIR=/glade/scratch/zarzycki/LES-snow/stats/
EXTRACTOUTDIR=/glade/scratch/zarzycki/LES-snow/storm-files/

PRECTFILE="/glade/p_old/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PRECT/b.e11."${RUNCONFIG}".f09_g16."${NUMID}".cam.h2.PRECT."${RUNDATES}".nc"
PRECVARNAME="PRECT"
PSLFILE="/glade/p_old/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PSL/b.e11."${RUNCONFIG}".f09_g16."${NUMID}".cam.h2.PSL."${RUNDATES}".nc"
PSLVARNAME="PSL"
SNOWFILE="/glade/u/home/zarzycki/scratch/LENS-snow/b.e11."${RUNCONFIG}".f09_g16."${NUMID}".cam.h2.PRECT_SNOW."${RUNDATES}".nc"
SNOWVARNAME="PRECT_SNOW"

TRAJFILE="./traj/conntraj."${DESCSTR}
EXTRACTOUTFILE="ind-storms."${DESCSTR}".nc"

DATESTRING=`date +"%s%N"`

TEMPESTFILE="./traj/tempest."${MODEL}"."${RUNDATES:0:4}"."${NUMID}
TMPTRAJNAME="./traj/trajectories.txt."${MODEL}"."${RUNDATES:0:4}"."${NUMID}

starttime=$(date -u +"%s")

###################################################################################################

if [ ${DO_TEMPEST} == 'true' ]; then

  starttimetrack=$(date -u +"%s")

  BUILDDIR=/glade/work/zarzycki/tempestextremes_noMPI/bin/

  rm cyc.${DATESTRING}
  rm ${TMPTRAJNAME}
  touch cyc.${DATESTRING}

  FILELISTNAME=filelist.txt.${DATESTRING}
  ls "/glade/p_old/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PSL/b.e11."${RUNCONFIG}".f09_g16."${NUMID}".cam.h2.PSL."${RUNDATES}".nc" -1 > ${FILELISTNAME}

  ${BUILDDIR}/DetectCyclonesUnstructured --in_data_list "${FILELISTNAME}" --out cyc_tempest.${DATESTRING} --closedcontourcmd "PSL,200.,6.,0" --outputcmd "PSL,min,0" --maxlat 80. --minlat 10. </dev/null
  cat cyc_tempest.${DATESTRING}* >> cyc.${DATESTRING}
  rm cyc_tempest.${DATESTRING}*

  rm ${FILELISTNAME}
  rm log*

  ${BUILDDIR}/StitchNodes --format "i,j,lon,lat,slp" --in cyc.${DATESTRING} --range 8.0 --minlength 5 --maxgap 0 --out ${TMPTRAJNAME} --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

  #BUILDDIR=/glade/p_old/work/zarzycki/tempestextremes/bin/
  #${BUILDDIR}/StitchNodes --format "i,j,lon,lat,slp" --in cyc.1505498447694611008 --range 8.0 --minlength 5 --maxgap 0 --out trajectories.txt.LENS.5e-9_12.001 --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

  endtimetrack=$(date -u +"%s")

  # Apparently need this to make sure stuff is written to disk?
  sleep 5

  starttimefilter=$(date -u +"%s")

  ncl filter_ne_storms.ncl 'infile="'${TMPTRAJNAME}'"' \
   'outfile="'${TEMPESTFILE}'"'

  rm cyc.${DATESTRING}
  rm ${TMPTRAJNAME}

  endtimefilter=$(date -u +"%s")

fi # end do tempest


if [ ${DO_CONNTRACK} == 'true' ]; then
  starttimefilter=$(date -u +"%s")

  ncl ${NCLPATH}/connected-tracking.ncl 'snowFileFull="'${PRECTFILE}'"' \
     'outfile="'${TRAJFILE}'"' \
     'thresh_prec_str="'${THRESH}'"'

  endtimetrack=$(date -u +"%s")
fi

###################################################################################################

if [ ${SEASONFRAC_TEMPEST} == 'true' ]; then
  OUTFRACFILE="frac.tempest."${DESCSTR}".nc"
  (set -x; ncl ${NCLPATH}/calc-seasonal-fraction.ncl 'snowFileFull="'${SNOWFILE}'"' \
    'SNOWVARNAME="'${SNOWVARNAME}'"' \
    'prectFileFull="'${PRECTFILE}'"' \
    'traj_filename="'${TEMPESTFILE}'"' \
    'outDir="'${EXTRACTOUTDIR}'"' \
    'outFileName="'${OUTFRACFILE}'"')
fi

if [ ${SEASONFRAC_CONNTRACK} == 'true' ]; then
  OUTFRACFILE="frac.conntraj."${DESCSTR}".nc"
  (set -x; ncl ${NCLPATH}/calc-seasonal-fraction.ncl 'snowFileFull="'${SNOWFILE}'"' \
    'SNOWVARNAME="'${SNOWVARNAME}'"' \
    'prectFileFull="'${PRECTFILE}'"' \
    'traj_filename="'${TRAJFILE}'"' \
    'outDir="'${EXTRACTOUTDIR}'"' \
    'outFileName="'${OUTFRACFILE}'"')
fi

###################################################################################################

if [ ${PROCESS_CONNTRACK} == 'true' ]; then

  starttimeextract=$(date -u +"%s")

  (set -x; ncl ${NCLPATH}/extract_individual_storms.ncl 'snowFileFull="'${SNOWFILE}'"' \
    'SNOWVARNAME="'${SNOWVARNAME}'"' \
    'prectFileFull="'${PRECTFILE}'"' \
    'PRECTVARNAME="'${PRECVARNAME}'"' \
    'pslFileFull="'${PSLFILE}'"' \
    'PSLVARNAME="'${PSLVARNAME}'"' \
    'traj_filename="'${TRAJFILE}'"' \
    'outDir="'${EXTRACTOUTDIR}'"' \
    'outFileName="'${EXTRACTOUTFILE}'"')

   endtimeextract=$(date -u +"%s")
   starttimeRSI=$(date -u +"%s")

  RSISNOWFILES=(${RSISNOWFILE}".SNOW.conntraj.csv" ${RSISNOWFILE}".PRECT.conntraj.csv")
  SWES=(${SWE} ${SWE})
  SNOWVARNAMES=("CUM_SNOWFALL" "CUM_PRECT")
  # LOOPVAR does both snow and total precip (pretending total precip is snow)
  for LOOPVAR in `seq 0 1`
  do
    rm ${RSIOUTDIR}/${RSISNOWFILES[${LOOPVAR}]} 
    LOOPST=0
    LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
    for IX in `seq ${LOOPST} ${LOOPEN}`
    do
      (set -x; ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} SWE=${SWES[${LOOPVAR}]} \
       'imgDir="'${RSIOUTDIR}'/images-RSI-tempest/'${MODEL}'.'${RUNDATES:0:4}'.'${NUMID}'"' \
       'RSIoutFile="'${RSIOUTDIR}/${RSISNOWFILES[${LOOPVAR}]}'"' \
       'stormFilePath="'${EXTRACTOUTDIR}/${EXTRACTOUTFILE}'"' \
       'configStr="'${DESCSTR}'"' \
       'SNOWVARNAME="'${SNOWVARNAMES[${LOOPVAR}]}'"')
    done
  done
  endtimeRSI=$(date -u +"%s")
fi

###################################################################################################

if [ ${PROCESS_TEMPEST} == 'true' ]; then
  starttimeextract=$(date -u +"%s")
  TRAJFILE=${TEMPESTFILE}
  EXTRACTOUTFILE=${EXTRACTOUTFILE}".tempest.nc"
  ncl ${NCLPATH}/extract_individual_storms.ncl 'snowFileFull="'${SNOWFILE}'"' \
    'SNOWVARNAME="'${SNOWVARNAME}'"' \
    'prectFileFull="'${PRECTFILE}'"' \
    'PRECTVARNAME="'${PRECVARNAME}'"' \
    'pslFileFull="'${PSLFILE}'"' \
    'PSLVARNAME="'${PSLVARNAME}'"' \
    'traj_filename="'${TRAJFILE}'"' \
    'outDir="'${EXTRACTOUTDIR}'"' \
    'outFileName="'${EXTRACTOUTFILE}'"'
  endtimeextract=$(date -u +"%s")

  starttimeRSI=$(date -u +"%s")
  RSISNOWFILES=(${RSISNOWFILE}".SNOW.tempest.csv" ${RSISNOWFILE}".PRECT.tempest.csv")
  SWES=(${SWE} ${SWE})
  SNOWVARNAMES=("CUM_SNOWFALL" "CUM_PRECT")
  # LOOPVAR does both snow and total precip (pretending total precip is snow)
  for LOOPVAR in `seq 0 1`
  do
    rm ${RSIOUTDIR}/${RSISNOWFILES[${LOOPVAR}]} 
    LOOPST=0
    LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
    for IX in `seq ${LOOPST} ${LOOPEN}`
    do
      (set -x; ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} SWE=${SWES[${LOOPVAR}]} \
       'imgDir="'${RSIOUTDIR}'/images-RSI-tempest/'${MODEL}'.'${RUNDATES:0:4}'.'${NUMID}'"' \
       'RSIoutFile="'${RSIOUTDIR}/${RSISNOWFILES[${LOOPVAR}]}'"' \
       'stormFilePath="'${EXTRACTOUTDIR}/${EXTRACTOUTFILE}'"' \
       'configStr="'${DESCSTR}'"' \
       'SNOWVARNAME="'${SNOWVARNAMES[${LOOPVAR}]}'"')
    done
  done
  endtimeRSI=$(date -u +"%s")
fi

# Calculate timings for each segment by checking if starttime was defined and then calculating if defined
if [ -z ${starttimetrack+x} ]; then tottimetrack=0; else tottimetrack=$(($endtimetrack-$starttimetrack)); fi
if [ -z ${starttimefilter+x} ]; then tottimefilter=0; else tottimefilter=$(($endtimefilter-$starttimefilter)); fi
if [ -z ${starttimeextract+x} ]; then tottimeextract=0; else tottimeextract=$(($endtimeextract-$starttimeextract)); fi
if [ -z ${starttimeRSI+x} ]; then tottimeRSI=0; else tottimeRSI=$(($endtimeRSI-$starttimeRSI)); fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

printf "${RUNDATES:0:4} ${NUMID} ${tottimetrack} ${tottimefilter} ${tottimeextract} ${tottimeRSI} ${tottime}\n" >> timing.txt
