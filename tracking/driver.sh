#!/bin/bash

##### USER OPTIONS
MACH=CORI     # current options are CORI, CHEYENNE, SERIAL

COMPONENTDIR="./gen-script/"
if [ "$MACH" == "CORI" ]; then
   echo "setting batch settings for CORI"
   MACHFILE=${COMPONENTDIR}/mach.cori
   SUBCMD="sbatch"
elif [ "$MACH" == "CHEYENNE" ]; then
   echo "setting batch settings for CHEYENNE"
   MACHFILE=${COMPONENTDIR}/mach.cheyenne
   SUBCMD="qsub"
elif [ "$MACH" == "SERIAL" ]; then
   echo "setting batch settings for SERIAL"
   MACHFILE=${COMPONENTDIR}/mach.login
   SUBCMD=""
else
   echo "Unknown machine, exiting..."
   exit
fi

echo $MACHFILE

cat $MACHFILE nl.SHELL.hyperion ${COMPONENTDIR}/mechanics.txt > TMP.sh

echo "submitting job"
${SUBCMD} TMP.sh
