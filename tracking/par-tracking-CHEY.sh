#!/bin/bash

##=======================================================================
#PBS -N tempest-ETCs
#PBS -A P54048000 
#PBS -l walltime=0:29:00
#PBS -q premium
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mpiprocs=36
################################################################

############ USER OPTIONS #####################

## Path to TempestExtremes binaries on YS
#TEMPESTEXTREMESDIR=~/software/tempestextremes/
TEMPESTEXTREMESDIR=/glade/work/zarzycki/tempestextremes/

#UQSTR=NE30
#TOPOFILE=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc
#PATHTOFILES=/glade2/h2/acgd0005/archive/f.asd2017.cesm20b05.FAMIPC6CLM5.ne30_g16.exp212/atm/hist/
#CONNECTFLAG="--in_connect /glade/u/home/zarzycki/tempest-scripts/asd-se/ne30.connect.dat"

UQSTR=MP15A-120A-US-EXP213
TOPOFILE=/glade/work/zarzycki/ASD2017_files/atm/cam/topo/mp15a-120a-US.topo.170118.nc
PATHTOFILES=/glade2/h2/acgd0005/archive/f.asd2017.cesm20b05.FAMIPC6CLM5.mp15a-120a-US_t12.exp213/atm/hist/
CONNECTFLAG="--in_connect /glade/u/home/zarzycki/tempest-scripts/asd-se/mp15a-120a-US.connect.dat"

#module load impi
MPIRUNCMD="mpiexec_mpt" # mpiexec_mpt   srun -n 32    mpirun -n 36

############ TRACKER MECHANICS #####################

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
TRAJFILENAME=trajectories.txt.${UQSTR}
touch $FILELISTNAME

# Generate parallel file
FILES=${PATHTOFILES}/*.h3.*.nc
for f in $FILES
do
  echo "${f};${TOPOFILE}" >> $FILELISTNAME
done

starttime=$(date -u +"%s")

STR_DETECT="--verbosity 0 --timestride 2 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --mergedist 8.0 --closedcontourcmd PSL,200.,6.,0 --outputcmd PSL,min,0"
touch cyclones.${DATESTRING}
${MPIRUNCMD} ${TEMPESTEXTREMESDIR}/bin/DetectCyclonesUnstructured --in_data_list "${FILELISTNAME}" ${STR_DETECT} </dev/null
cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
rm cyclones_tempest.${DATESTRING}*

# Stitch candidate cyclones together
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "ncol,lon,lat,slp" --in cyclones.${DATESTRING} --range 8.0 --minlength 5 --maxgap 0 --out ${TRAJFILENAME} --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

printf "${tottime}\n" >> timing.txt

rm cyclones.${DATESTRING}
rm ${FILELISTNAME}
rm log*
