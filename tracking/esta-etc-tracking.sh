#!/bin/bash

##=======================================================================
#PBS -N tempest-ETCs
#PBS -A P05010048 
#PBS -l walltime=0:29:00
#PBS -q premium
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mpiprocs=36
################################################################

#qsub -I -l select=1:ncpus=36:mpiprocs=36 -l walltime=01:00:00 -q premium -A P05010048

if [ -z "$1" ]; then
  echo "No namelist passed to script"
  exit
elif [ ! -f "$1" ]; then
  echo "$1 does not exist."
  exit
else
  echo "Using $1"
fi

NAMELISTFILE=${1}
while IFS= read -r line; do
  IN="$line"
  arrIN=(${IN//=/ })
  VARNAME=${arrIN[0]}
  VALUE=`echo "${arrIN[1]}" | cut -f1 -d","`
  echo "NAMELIST: setting ${VARNAME} to ${VALUE}"
  eval $VARNAME=$VALUE
done < "$NAMELISTFILE"

export ESTAPATH  # need to export such that filter_ne_storms can see

############ TRACKER MECHANICS #####################

DATESTRING=`date +"%s%N"`
TRAJFILENAME=tmp.trajectories.txt.${TE_UQSTR}
TRAJFILE=${TRAJDIR}"/traj."${TE_UQSTR}

starttime=$(date -u +"%s")

STR_DETECT="--verbosity 0 --timestride 1 ${TE_CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --mergedist 8.0 --closedcontourcmd PSL,200.,6.,0 --outputcmd PSL,min,0"
touch cyclones.${DATESTRING}
${TE_MPIRUNCMD} ${TE_TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STR_DETECT} </dev/null
cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
rm cyclones_tempest.${DATESTRING}*

# Stitch candidate cyclones together
${TE_TEMPESTEXTREMESDIR}/bin/StitchNodes --format "i,j,lon,lat,slp" --in cyclones.${DATESTRING} --range 8.0 --minlength 5 --maxgap 0 --out ${TRAJFILENAME} --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

ncl ${ESTAPATH}/tracking/filter_ne_storms.ncl 'infile="'${TRAJFILENAME}'"' 'outfile="'${TRAJFILE}'"'

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

printf "`basename $0` ${TE_UQSTR} ${tottime}\n" >> timing.txt

# Cleanup
rm cyclones.${DATESTRING}
rm log*
rm ${TRAJFILENAME}
