#!/bin/bash

##=======================================================================
#PBS -N gnu-par
#PBS -A P54048000 
#PBS -l walltime=11:59:00
#PBS -q regular
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mem=109GB
################################################################

date

module load parallel

NUMCORES=18
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt
rm ${COMMANDFILE}

declare -a RUNDATES_S=("1990010100Z-2005123118Z" "2026010100Z-2035123118Z" "2071010100Z-2080123118Z")
#declare -a RUNDATES_S=("1990010100Z-2005123118Z")
declare -a LWE_S=("13")
declare -a DELTALAT_S=("15")
# create new line in command file for 1 -> NUMCORES
#for i in `seq 1 ${NUMTASKS}`;

for LWE in "${LWE_S[@]}"
do
  for RUNDATES in "${RUNDATES_S[@]}"
  do
    for DELTALAT in "${DELTALAT_S[@]}"
    do
      for ii in {001..035};       # 001..042 processes all available LENS data
      do
        LINECOMMAND="./lens-snow-tracker.sh ${ii} \"${RUNDATES}\" ${LWE} 5e-9 ${DELTALAT}  " 
        echo ${LINECOMMAND} >> ${COMMANDFILE}
      done
    done
  done
done

# Launch GNU parallel
#### Use this for login nodes (nohup ./batch.sh &)
#parallel --jobs ${NUMCORES} < ${COMMANDFILE}

#### Use this for Cheyenne batch jobs
parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}

date 
