#!/bin/bash

##=======================================================================
#PBS -N gnu-par
#PBS -A P54048000 
#PBS -l walltime=11:00:00
#PBS -q regular
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36
################################################################

date

module load parallel

NUMCORES=4
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

rm ${COMMANDFILE}
# create new line in command file for 1 -> NUMCORES
#for i in `seq 1 ${NUMTASKS}`;
for ii in {1958..2017};
do
  LINECOMMAND="./jra-snow-tracker.sh ${ii}  " 
  echo ${LINECOMMAND} >> ${COMMANDFILE}
done

# Launch GNU parallel
parallel --jobs ${NUMCORES} < ${COMMANDFILE}

rm ${COMMANDFILE}

date 
