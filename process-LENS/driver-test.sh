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
ENSNUM=""
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
