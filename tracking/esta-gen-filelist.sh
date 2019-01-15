#!/bin/bash

#########################

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

########################

DATESTRING=`date +"%s%N"`
rm $FILELISTNAME ; touch $FILELISTNAME
# Generate parallel file
FILES=${TE_PATHTOFILES}
for f in $FILES
do
  echo $f
  echo "${f}" >> $FILELISTNAME
done
