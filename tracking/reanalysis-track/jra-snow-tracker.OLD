#!/bin/bash

##=======================================================================
#BSUB -a poe                     # use LSF openmp elim
#BSUB -N
#BSUB -n 1                      # yellowstone setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q geyser                 # queue
#BSUB -J sub_ncl 
#BSUB -W 6:00                   # wall clock limit
#BSUB -P P05010048               # account number

################################################################

date

  NUMID=$1

  DO_TEMPEST=true
  DO_CONNTRACK=false
  PROCESS_TEMPEST=true
  PROCESS_CONNTRACK=false

  SWE="12"
  THRESH="5e-9"

#b.e11.B20TRC5CNBDRD.f09_g16.027.cam.h2.PRECT_SNOW_RATE.1990010100Z-2005123118Z.n
  MODEL="JRA"
  RUNDATES=${NUMID}
  RUNCONFIG="X"
  CONFIG=${THRESH}"_"${SWE}

  DESCSTR=${MODEL}"."${RUNDATES:0:4}"."${NUMID}"."${CONFIG}

  RSIDIR=/glade/u/home/zarzycki/snow-tracking/calc_RSI/
  RSISNOWFILE="RSI.SNOW."${DESCSTR}".csv"
  RSIPRECFILE="RSI.PREC."${DESCSTR}".csv"

  NCLPATH=/glade/u/home/zarzycki/snow-tracking/tracking/
  RSIOUTDIR=/glade/scratch/zarzycki/LES-snow/stats/
  EXTRACTOUTDIR=/glade/scratch/zarzycki/LES-snow/storm-files/

  PRECTFILE="/glade/u/home/zarzycki/scratch/h1files/JRA/"${NUMID}"/JRA.h1."${NUMID}".PRECT.nc"
  PRECVARNAME="PRECT"
  PSLFILE="/glade/u/home/zarzycki/scratch/h1files/JRA/"${NUMID}"/JRA.h1."${NUMID}".PSL.nc"
  PSLVARNAME="PSL"
  SNOWFILE="/glade/u/home/zarzycki/scratch/h1files/JRA/"${NUMID}"/JRA.h1."${NUMID}".PRECSN.nc"
  SNOWVARNAME="PRECSN"

  TRAJFILE="./traj/conntraj."${DESCSTR}
  EXTRACTOUTFILE="ind-storms."${DESCSTR}".nc"

  DATESTRING=`date +"%s%N"`

  TEMPESTFILE="./traj/tempest."${MODEL}"."${RUNDATES:0:4}"."${NUMID}
  TMPTRAJNAME="./traj/trajectories.txt."${MODEL}"."${RUNDATES:0:4}"."${NUMID}

###################################################################################################

  if [ ${DO_TEMPEST} == 'true' ]; then
    BUILDDIR=/glade/p/work/zarzycki/tempestextremes/bin/

    rm cyc.${DATESTRING}
    rm ${TMPTRAJNAME}
    touch cyc.${DATESTRING}

    FILELISTNAME=filelist.txt.${DATESTRING}
    ls ${PSLFILE} -1 > ${FILELISTNAME}

    ${BUILDDIR}/DetectCyclonesUnstructured --in_data_list "${FILELISTNAME}" --out cyc_tempest.${DATESTRING} --closedcontourcmd "PSL,200.,6.,0" --outputcmd "PSL,min,0" --mergedist 8.0 --maxlat 80. --minlat 10. </dev/null
    cat cyc_tempest.${DATESTRING}* >> cyc.${DATESTRING}
    rm cyc_tempest.${DATESTRING}*

    rm ${FILELISTNAME}
    rm log*

    ${BUILDDIR}/StitchNodes --format "i,j,lon,lat,slp" --in cyc.${DATESTRING} --range 8.0 --minlength 5 --maxgap 0 --out ${TMPTRAJNAME} --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

    #BUILDDIR=/glade/p/work/zarzycki/tempestextremes/bin/
    #${BUILDDIR}/StitchNodes --format "i,j,lon,lat,slp" --in cyc.1505498447694611008 --range 8.0 --minlength 5 --maxgap 0 --out trajectories.txt.LENS.5e-9_12.001 --threshold "lat,>,25,3;lon,>,260,3;lat,<,65,3;lon,<,310,3"

    endtimetrack=$(date -u +"%s")

    # Apparently need this to make sure stuff is written to disk?
    sleep 5

    ncl filter_ne_storms.NEW.ncl 'infile="'${TMPTRAJNAME}'"' \
     'outfile="'${TEMPESTFILE}'"'

    rm cyc.${DATESTRING}
    rm ${TMPTRAJNAME}

  fi # end do tempest


  if [ ${DO_CONNTRACK} == 'true' ]; then
    ncl ${NCLPATH}/connected-tracking.ncl 'snowFileFull="'${PRECTFILE}'"' \
       'outfile="'${TRAJFILE}'"' \
       'thresh_prec_str="'${THRESH}'"'
  fi

#############

  if [ ${PROCESS_CONNTRACK} == 'true' ]; then

    ncl ${NCLPATH}/extract_individual_storms.ncl 'snowFileFull="'${SNOWFILE}'"' \
      'SNOWVARNAME="'${SNOWVARNAME}'"' \
      'prectFileFull="'${PRECTFILE}'"' \
      'PRECTVARNAME="'${PRECVARNAME}'"' \
      'pslFileFull="'${PSLFILE}'"' \
      'PSLVARNAME="'${PSLVARNAME}'"' \
      'traj_filename="'${TRAJFILE}'"' \
      'outDir="'${EXTRACTOUTDIR}'"' \
      'outFileName="'${EXTRACTOUTFILE}'"' \
      'SWE_str="'${SWE}'"'
 
     rm ${RSIOUTDIR}/${RSISNOWFILE}
     LOOPST=0
     LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
     for IX in `seq ${LOOPST} ${LOOPEN}`
     do
       ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} \
        'imgDir="'${RSIOUTDIR}'/images-RSI-conntraj/"' \
        'RSIoutFile="'${RSIOUTDIR}/${RSISNOWFILE}'"' \
        'stormFilePath="'${EXTRACTOUTDIR}/${EXTRACTOUTFILE}'"' \
        'configStr="'${DESCSTR}'"' \
        'SNOWVARNAME="CUM_SNOWFALL"'
     done

  fi

####################

  if [ ${PROCESS_TEMPEST} == 'true' ]; then

    TRAJFILE=${TEMPESTFILE}
    EXTRACTOUTFILE=${EXTRACTOUTFILE}".tempest.nc"

    ncl ${NCLPATH}/extract_individual_storms.ncl 'snowFileFull="'${SNOWFILE}'"' \
      'SNOWVARNAME="'${SNOWVARNAME}'"' \
      'prectFileFull="'${PRECTFILE}'"' \
      'PRECTVARNAME="'${PRECVARNAME}'"' \
      'pslFileFull="'${PSLFILE}'"' \
      'PSLVARNAME="'${PSLVARNAME}'"' \
      'traj_filename="'${TRAJFILE}'"' \
      'outDir="'${EXTRACTOUTDIR}'"' \
      'outFileName="'${EXTRACTOUTFILE}'"'

    RSISNOWFILES=(${RSISNOWFILE}".SNOW.tempest.csv" ${RSISNOWFILE}".PRECT.tempest.csv")
    SWES=(${SWE} ${SWE})
    SNOWVARNAMES=("CUM_SNOWFALL" "CUM_PRECT")

    for LOOPVAR in `seq 0 1`
    do
      rm ${RSIOUTDIR}/${RSISNOWFILE}
      LOOPST=0
      LOOPEN=$((`grep -c start ${TRAJFILE}` - 1))
      for IX in `seq ${LOOPST} ${LOOPEN}`
      do
        ncl ${RSIDIR}/calculate_RSI.ncl stormID=${IX} SWE=${SWES[${LOOPVAR}]} \
         'imgDir="'${RSIOUTDIR}'/images-RSI-tempest/"' \
         'RSIoutFile="'${RSIOUTDIR}/${RSISNOWFILES[${LOOPVAR}]}'"' \
         'stormFilePath="'${EXTRACTOUTDIR}/${EXTRACTOUTFILE}'"' \
         'configStr="'${DESCSTR}'"' \
         'SNOWVARNAME="'${SNOWVARNAMES[${LOOPVAR}]}'"'
      done
    done
  fi

date 
