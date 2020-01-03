#!/bin/bash

NAMELISTFILE=${1}

while IFS= read -r line; do
  IN="$line"
  arrIN=(${IN//=/ })
  VARNAME=${arrIN[0]}
  VALUE=`echo "${arrIN[1]}" | cut -f1 -d","`
  echo "NAMELIST: setting ${VARNAME} to ${VALUE}"
  eval $VARNAME=$VALUE
done < "$NAMELISTFILE"