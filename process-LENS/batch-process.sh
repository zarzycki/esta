#!/bin/bash

##=======================================================================
#PBS -N gnu-par
#PBS -A P54048000 
#PBS -l walltime=06:00:00
#PBS -q regular
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mem=109GB
################################################################

NUMCORES=15
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

CAMDIR=/glade/scratch/zarzycki/ASD/archive/f.asd2017.cesm20b05.FAMIPC6CLM5.mp120a_g16.exp214/atm/hist/
#CAMDIR=/glade/u/home/zarzycki/acgd0005/archive/f.asd2017.cesm20b05.FAMIPC6CLM5.ne30_g16.exp212/atm/hist/
OUTDIR=/glade/scratch/zarzycki/TEST-SNOW/
FILES=`ls ${CAMDIR}/*.cam.h2.*.nc`
for f in $FILES
  do
  echo "Processing $f file..."
  base=`basename $f`
  OUTFILE=${base/h2/h5}
  LINECOMMAND="./1out ${f} ${OUTDIR}/${OUTFILE} 0"
  echo ${LINECOMMAND} >> ${COMMANDFILE}
  OUTFILE=${base/h2/h6}
  LINECOMMAND="./1out ${f} ${OUTDIR}/${OUTFILE} 1"
  echo ${LINECOMMAND} >> ${COMMANDFILE}
  OUTFILE=${base/h2/h7}
  LINECOMMAND="./1out ${f} ${OUTDIR}/${OUTFILE} 2"
  echo ${LINECOMMAND} >> ${COMMANDFILE}
  OUTFILE=${base/h2/h8}
  LINECOMMAND="./1out ${f} ${OUTDIR}/${OUTFILE} 3"
  echo ${LINECOMMAND} >> ${COMMANDFILE}
  #./1out ${f} ${OUTDIR}/${OUTFILE}
  OUTFILE=${base/h2/h2}
  LINECOMMAND="ncks -v PRECBSN,PRECBRA,PRECBIP,PRECBFZ,PRECT $f ${OUTDIR}/${OUTFILE}"
  echo ${LINECOMMAND} >> ${COMMANDFILE}
done

#### Use this for Cheyenne batch jobs
parallel --jobs ${NUMCORES} -u --sshloginfile $PBS_NODEFILE --workdir $PWD < ${COMMANDFILE}

rm ${COMMANDFILE}
