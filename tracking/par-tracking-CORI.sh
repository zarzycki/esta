#!/bin/bash -l

################################################################
#SBATCH -N 1                #Number of nodes
#SBATCH -t 00:30:00         #Time limit
#SBATCH -q premium          #Use the regular QOS
#SBATCH -L SCRATCH          #Job requires $SCRATCH file system
#SBATCH -C knl,quad,cache   #Use KNL nodes in quad cache format (default, recommended)
################################################################

# salloc -N 1 -C knl -q interactive -t 01:00:00

############ USER OPTIONS #####################

## Path to TempestExtremes binaries on YS
TEMPESTEXTREMESDIR=~/software/tempestextremes/

CONFIG=dtime900.003
MACH=CORI
GRID=WAT
UQSTR=${MACH}.VR28.NATL.${GRID}.CAM5.4CLM5.0.${CONFIG}
TOPOFILE=/global/homes/c/czarzyck/scratch/unigridFiles/ne0np4natlantic${GRID,,}.ne30x4/topo/topo_ne0np4natlantic${GRID,,}.ne30x4_interpic.nc
PATHTOFILES=/global/homes/c/czarzyck/scratch/hyperion/${MACH}.VR28.NATL.${GRID}.CAM5.4CLM5.0.${CONFIG}/atm/hist/
CONNECTFLAG="--in_connect /global/homes/c/czarzyck/tempest-scripts/hyperion/ne0np4natlantic${GRID,,}.ne30x4.connect.dat" 

############ TRACKER MECHANICS #####################

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
TRAJFILENAME=trajectories.txt.${UQSTR}
FILTFILE=traj.filt
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
srun -n 32 ${TEMPESTEXTREMESDIR}/bin/DetectCyclonesUnstructured --in_data_list "${FILELISTNAME}" ${STR_DETECT} </dev/null
cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
rm cyclones_tempest.${DATESTRING}*

# Stitch candidate cyclones together
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "ncol,lon,lat,slp" --in cyclones.${DATESTRING} --range 8.0 --minlength 5 --maxgap 0 --out ${TRAJFILENAME} --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

ncl filter_ne_storms.ncl 'infile="'${TRAJFILENAME}'"' 'outfile="'${FILTFILE}'"'

rm cyclones.${DATESTRING}
rm ${FILELISTNAME}
rm log*

