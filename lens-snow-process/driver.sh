#!/bin/bash

##=======================================================================
#BSUB -a poe                     # use LSF openmp elim
#BSUB -N
#BSUB -n 1                      # yellowstone setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q geyser                 # queue
#BSUB -J process_lens
#BSUB -W 11:59                   # wall clock limit
#BSUB -P P54048000               # account number

################################################################

starttime=$(date -u +"%s")
#ENSNUM="011"
TIMEPERIOD="1990010100Z-2005123118Z"  # 1990010100Z-2005123118Z, 2026010100Z-2035123118Z, 2071010100Z-2080123118Z
LENSCONFIG="B20TRC5CNBDRD"            # B20TRC5CNBDRD, BRCP85C5CNBDRD
INDIR="/glade/scratch/zarzycki/LES_snow/"
OUTDIR="/glade/scratch/zarzycki/LES_snow/"
ARCHIVEDIR="/glade/scratch/zarzycki/LENS-snow-proc/"
SCRIPTDIR=$PWD
declare -a ALLVARS=("Q" "T" "U" "V" "Z3" "PRECT" "PS")

# If ENSNUM is unset (or empty string), draw from members.txt
if ([ -z "$ENSNUM" ] && [ "${ENSNUM+xxx}" = "xxx" ]) || [ -z "${ENSNUM+xxx}" ] ; then
  echo "ENSNUM pulled from members.txt file"
  ENSNUM=$(head -n 1 members.txt)
  #Remove top line from members file
  tail -n +2 members.txt > members2.txt
  mv -v members2.txt members.txt
else
  echo "ENSNUM specified by user"
fi
echo "ENSNUM: "$ENSNUM

###Get HSI files
cd ${INDIR}
for ii in "${ALLVARS[@]}"
do
  VAR=${ii}
  if [ ! -f b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.${VAR}.${TIMEPERIOD}.nc ] ; then
    hsi cget /CCSM/csm/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/${VAR}/b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.${VAR}.${TIMEPERIOD}.nc
  fi
done

# Get number of time indices, decide how many files we need to split into
TOTLINES=`ncks -C -v time -m b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.${VAR}.${TIMEPERIOD}.nc | grep "size =" | grep "time, " | cut -d' ' -f 7`
LINEPERFILE=40
let NUMLOOP=($TOTLINES+$LINEPERFILE-1)/$LINEPERFILE; echo $NUMLOOP
#NUMLOOP=2

###Split LENS 6-hourly files for memory purposes
cd ${INDIR}
for ii in "${ALLVARS[@]}"
do
  VAR=${ii}
  for i in $(seq 1 $NUMLOOP)
  do
    START=$(( (i-1)*LINEPERFILE ))
    END=$(( (i*LINEPERFILE)-1 ))
    if [ "$i" -eq "$NUMLOOP" ] ; then
      END=$(( TOTLINES-1 ))
    fi
    printf -v INTSTR "%06d" $START
    echo ${ENSNUM} .. ${VAR} .. $START .. $END --- $INTSTR
    ncks -d time,$START,$END b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.${VAR}.${TIMEPERIOD}.nc ${VAR}_${ENSNUM}_${TIMEPERIOD}_${INTSTR}.nc
  done
done

###Process files
cd $SCRIPTDIR
for i in $(seq 1 $NUMLOOP)
#for i in $(seq 1 1)
do
  START=$(( (i-1)*LINEPERFILE )) ; printf -v INTSTR "%06d" $START
  ./1out ${ENSNUM} ${INTSTR} ${TIMEPERIOD} ${INDIR} ${OUTDIR}
done

cd ${OUTDIR}

# Remove other variables split
for ii in "${ALLVARS[@]}"
do
  ls ${ii}_${ENSNUM}_*.nc
  rm ${ii}_${ENSNUM}_*.nc
done

# Concat all SNOW files into one big time series
ncrcat SNOW_${ENSNUM}_*.nc b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.WINTER.${TIMEPERIOD}.nc
rm SNOW_${ENSNUM}_*.nc

# Extract variables from big time series
declare -a EXVARS=("PRECT_SNOW_RATE" "PRECT_FZRA" "PRECT_ICE" "PRECT_SNOW" "PRECT_RAIN" "PTYPE" "PTYPECZ" "RATIO")
for zz in "${EXVARS[@]}"
do
  EXVAR=${zz}
  echo "Extracting "${EXVAR}"..."
  ncks -v ${EXVAR} b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.WINTER.${TIMEPERIOD}.nc b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.${EXVAR}.${TIMEPERIOD}.nc
done

rm b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.WINTER.${TIMEPERIOD}.nc

for xx in "${ALLVARS[@]}"
do
  VAR=${xx}
  rm b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.${VAR}.${TIMEPERIOD}.nc
done

mkdir -p ${ARCHIVEDIR}
mv b.e11.${LENSCONFIG}.f09_g16.${ENSNUM}.cam.h2.*.${TIMEPERIOD}.nc ${ARCHIVEDIR}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
echo $tottime
