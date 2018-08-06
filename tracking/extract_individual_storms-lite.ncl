;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/global/homes/c/czarzyck/ncl-zarzycki/functions/getTrajectories.ncl"
load "fcns-extract-storm.ncl"
;========================================================

begin

wcStrt     = systemfunc("date")

;----------------Constants-----------------------------------

ms_to_kts = 1.94384449
pi = 3.14159265358979
omega = 7.2921E-5
ms_to_inhr=141732.

; NL VARS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;list_of_files="/global/homes/c/czarzyck/scratch/hyperion/CORI.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.003/atm/hist/*cam.h8.1986*.PTYPES.nc_remap.nc"
;VARNAME="PRECB_SN"
;traj_filename="./traj"
;outDir="/global/homes/c/czarzyck/scratch/hyperion/"
;outFileName="TEST.nc"
;AGGRO_OPERATOR="sum"  ; sum, avg, min, max
;offset=6*ms_to_inhr
;deltaDeg=15.0
;doTimeSeries=True
;doRegionalOutput=False
;forceDeleteOutput=True

list_of_files="/global/homes/c/czarzyck/scratch/hyperion/CORI.VR28.NATL.WAT.CAM5.4CLM5.0.dtime900.003/atm/hist/*cam.h3.1986*_remap.nc"
VARNAME="FLUT"
traj_filename="./traj"
outDir="/global/homes/c/czarzyck/scratch/hyperion/"
outFileName="TEST.nc"
AGGRO_OPERATOR="min"  ; sum, avg, min, max
offset=1.
deltaDeg=200.0
doTimeSeries=True
doRegionalOutput=False
forceDeleteOutput=False
; NL VARS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


;----------------Find min/max bounds-----------------------------------

if (doRegionalOutput) then   ; If yes, truncation output
  minLatOut=20.
  maxLatOut=55.
  minLonOut=-95.
  maxLonOut=-60.
else   ; if not, just set to very large neg/pos values since NCL indexes min/max auto this way
  minLatOut=-99999.
  maxLatOut=99999.
  minLonOut=-99999.
  maxLonOut=99999.
end if

;----------------Pack min/max lat/lon into single 4x1 var-----------------------------------

latBounds=packLatLonBounds(minLatOut,maxLatOut,minLonOut,maxLonOut)

;-------------Get trajectories------------------------------

thefile=traj_filename
ncols = 8                     ; assign number of columns
maxLines = -1
traj_data = getTrajectories(thefile,ncols,maxLines)

traj_lon = traj_data(:,1,:)
traj_lat = traj_data(:,2,:)
traj_pres = traj_data(:,3,:)/100.
traj_yyyy = traj_data(:,4,:)
traj_mm = traj_data(:,5,:)
traj_dd = traj_data(:,6,:) 
traj_hh = traj_data(:,7,:)

ndims = dimsizes(traj_lon)
nstorms = ndims(0)
storm_index=ispan(0,nstorms-1,1)
ntimes = ndims(1)

;-------------lookup table------------------------------

; Here we build a "lookup" table to find what times are on what files...
; get list of h3 files
filelist = systemfunc ("ls "+list_of_files)
timeArr=buildTimeLookup(filelist)

;-------------get relevant coords------------------------------
retlist=getCoords(filelist,latBounds)
lat=retlist[0]
lon=retlist[1]
timeunits=retlist[2]
opt=retlist[3]
delete(retlist)

;-------------setup arrays------------------------------

wcStrtExtract     = systemfunc("date")

retlist=setupArrays(nstorms,ntimes,lat,lon)
lat2D=retlist[0]
lon2D=retlist[1]
SNOWFALL=retlist[2]
sttime=retlist[3]
entime=retlist[4]
durtime=retlist[5]
delete(retlist)

;-------------extract data from files------------------------------

SNOWFALL=loopAndExtract(SNOWFALL,VARNAME,lat2D,lon2D,deltaDeg,nstorms,ntimes,sttime,entime,durtime,traj_lon,traj_lat,traj_yyyy,traj_mm,traj_dd,traj_hh,timeArr,timeunits,opt,latBounds,offset)

;-------------do summation------------------------------

if (AGGRO_OPERATOR .eq. "sum") then
  CUM_SNOWFALL=dim_sum_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="SUM_"+VARNAME
else if (AGGRO_OPERATOR .eq. "avg") then
  CUM_SNOWFALL=dim_avg_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="AVG_"+VARNAME
else if (AGGRO_OPERATOR .eq. "min") then
  CUM_SNOWFALL=dim_min_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="MIN_"+VARNAME
else if (AGGRO_OPERATOR .eq. "max") then
  CUM_SNOWFALL=dim_max_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="MAX_"+VARNAME
else
  print("invalid aggregator operator, exiting...")
  exit
end if
end if
end if
end if

wallClockElapseTime(wcStrtExtract, "Time to extract all snowfall events", 0)

;-------------write output file------------------------------

wcStrtWriteout     = systemfunc("date")

setfileoption("nc","Format","LargeFile")

outFile=outDir+"/"+outFileName

print("writing:   "+outFile)

system("mkdir -p "+outDir)
if (forceDeleteOutput) then
  print("Force deleting any existing file...")
  system("/bin/rm -f "+outFile)
  fileMode="c"
else ; check if file exists, if does, add as read, if not, create
  if (fileexists(outFile)) then
    print("File exists, appending...")
    fileMode="w"
  else
    print("File doesn't exist, creating...")
    fileMode="c"
  end if
end if

ncdf = addfile(outFile,fileMode)

if (fileMode .eq. "c")
  fAtt               = True            ; assign file attributes
  fAtt@Conventions   = "None"
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes
end if

; print 1-D vars
ncdf->sttime=sttime
ncdf->entime=entime
ncdf->durtime=durtime

; print 2-D vars
ncdf->$AGGRO_OUT_NAME$ = CUM_SNOWFALL(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})

; print 3-D vars
if (doTimeSeries) then
  ncdf->$VARNAME$ = SNOWFALL(:,:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
end if

wallClockElapseTime(wcStrtWriteout, "Time to write output", 0)
wallClockElapseTime(wcStrt, "Total run time", 0)













end
