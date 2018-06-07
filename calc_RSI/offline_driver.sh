#!/bin/bash

ENSMEM=XXX

LOOPST=0
LOOPEN=399

date1=$(date -u +"%s")
for IX in `seq ${LOOPST} ${LOOPEN}`
do
  ncl calculate_RSI.ncl stormID=${IX} 'ensmember="'${ENSMEM}'"'
done
date2=$(date -u +"%s")
diff=$(($date2-$date1))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."



