#!/bin/bash

##=======================================================================
#PBS -N gnu-par
#PBS -A P54048000 
#PBS -l walltime=02:00:00
#PBS -q regular
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mem=109GB
################################################################

date

module load parallel

NUMCORES=36
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

rm ${COMMANDFILE}
# create new line in command file for 1 -> NUMCORES
#for i in `seq 1 ${NUMTASKS}`;
for ii in {001..035};       # 001..042 processes all available LENS data
do
  LINECOMMAND="./lens-snow-tracker.sh ${ii}  " 
  echo ${LINECOMMAND} >> ${COMMANDFILE}
done

# Launch GNU parallel
#### Use this for login nodes (nohup ./batch.sh &)
#parallel --jobs ${NUMCORES} < ${COMMANDFILE}

#### Use this for Cheyenne batch jobs
parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}

date 
