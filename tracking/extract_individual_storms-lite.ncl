;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

ESTAPATH = getenv("ESTAPATH")
fcntoload=ESTAPATH+"/functions/getTrajectoriesGeneric.ncl"
if (fileexists(fcntoload)) then
  loadscript(fcntoload)
else
  print("function "+fcntoload+" doesn't exist")
  exit
end if
fcntoload=ESTAPATH+"/functions/fcns-extract-storm.ncl"
if (fileexists(fcntoload)) then
  loadscript(fcntoload)
else
  print("function "+fcntoload+" doesn't exist")
  exit
end if
fcntoload=ESTAPATH+"/functions/namelist_extract.ncl"
if (fileexists(fcntoload)) then
  loadscript(fcntoload)
else
  print("function "+fcntoload+" doesn't exist")
  exit
end if

;========================================================

begin

wcStrt     = systemfunc("date")

;----------------Constants-----------------------------------

ms_to_kts = 1.94384449
pi = 3.14159265358979
omega = 7.2921E-5
ms_to_inhr=141732.

;nlfile="nl.gefs"

traj_dir = namelist_extract( nlfile, "TRAJDIR", "string")
UQSTR = namelist_extract( nlfile, "TE_UQSTR", "string")
outDir = namelist_extract( nlfile, "EX_OUTDIR", "string")
outFileName = namelist_extract( nlfile, "EX_OUTFILE", "string")
deltaDeg = namelist_extract( nlfile, "DELTADEG", "float")
list_of_files = namelist_extract( nlfile, "LIST_SNOW", "string")
VARNAME = namelist_extract( nlfile, "EX_INVARS", "string")
AGGRO_OPERATOR = namelist_extract( nlfile, "EX_AGGOP", "string")
offset = namelist_extract( nlfile, "EX_OFFSET", "float")
doTimeSeries = namelist_extract( nlfile, "EX_DOTIMESERIES", "logical")
doRegionalOutput = namelist_extract( nlfile, "EX_DOREGOUT", "logical")
forceDeleteOutput = namelist_extract( nlfile, "EX_FORCEDELOUT", "logical")

;-------------Get trajectories------------------------------

thefile=traj_dir+"/traj."+UQSTR
ncols = 9                     ; assign number of columns
maxLines = -1
headerDelimStr = "start"
delimStr="tab"
isUnstruc=False
isHeader=False
traj_data = getTrajectoriesGeneric(thefile,ncols,maxLines,headerDelimStr,delimStr,isUnstruc,isHeader)

traj_lon = traj_data(:,2,:)
traj_lat = traj_data(:,3,:)
traj_pres = traj_data(:,4,:)/100.
traj_yyyy = traj_data(:,5,:)
traj_mm = traj_data(:,6,:)
traj_dd = traj_data(:,7,:) 
traj_hh = traj_data(:,8,:)

delete(traj_data)

ndims = dimsizes(traj_lon)
nstorms = ndims(0)
storm_index=ispan(0,nstorms-1,1)
ntimes = ndims(1)

do ii=0,dimsizes(VARNAME)-1

;----------------Find min/max bounds-----------------------------------

if (doRegionalOutput(ii)) then   ; If yes, truncation output
  ;minLatOut=25.
  ;maxLatOut=55.
  ;minLonOut=265.
  ;maxLonOut=300.
  minLatOut=25.
  maxLatOut=55.
  minLonOut=265.;-360
  maxLonOut=300.;-360
else   ; if not, just set to very large neg/pos values since NCL indexes min/max auto this way
  minLatOut=-99999.
  maxLatOut=99999.
  minLonOut=-99999.
  maxLonOut=99999.
end if

;----------------Pack min/max lat/lon into single 4x1 var-----------------------------------

latBounds=packLatLonBounds(minLatOut,maxLatOut,minLonOut,maxLonOut)

;-------------lookup table------------------------------

; Here we build a "lookup" table to find what times are on what files...
; get list of h3 files
filelist = asciiread(list_of_files(ii),-1,"string")
;filelist = systemfunc("ls "+list_of_files(ii))  
timeArr=buildTimeLookup(filelist)

;-------------get relevant coords------------------------------
retlist=getCoords(filelist,latBounds)
lat=tofloat(retlist[0])
lon=tofloat(retlist[1])
if (minLonOut .ne. retlist[2]) then
  minLonOut=retlist[2]
  maxLonOut=retlist[3]
  print("repacking...")
  latBounds=packLatLonBounds(minLatOut,maxLatOut,minLonOut,maxLonOut)
end if
timeunits=retlist[4]
opt=retlist[5]
delete(retlist)
delete(filelist)

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

;list_vars()

;-------------extract data from files------------------------------

SNOWFALL=loopAndExtract(SNOWFALL,VARNAME(ii),lat2D,lon2D,deltaDeg,nstorms,ntimes,sttime,entime,durtime,traj_lon,traj_lat,traj_yyyy,traj_mm,traj_dd,traj_hh,timeArr,timeunits,opt,latBounds,offset(ii))

delete(timeArr)

;-------------do summation------------------------------

if (AGGRO_OPERATOR(ii) .eq. "sum") then
  CUM_SNOWFALL=dim_sum_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="SUM_"+VARNAME(ii)
else if (AGGRO_OPERATOR(ii) .eq. "avg") then
  CUM_SNOWFALL=dim_avg_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="AVG_"+VARNAME(ii)
else if (AGGRO_OPERATOR(ii) .eq. "min") then
  CUM_SNOWFALL=dim_min_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="MIN_"+VARNAME(ii)
else if (AGGRO_OPERATOR(ii) .eq. "max") then
  CUM_SNOWFALL=dim_max_n_Wrap(SNOWFALL,1)
  AGGRO_OUT_NAME="MAX_"+VARNAME(ii)
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

setfileoption("nc","Format","NetCDF4") 

outFile=outDir+"/"+outFileName

print("writing:   "+outFile)

system("mkdir -p "+outDir)
if (forceDeleteOutput(ii)) then
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
if (doTimeSeries(ii)) then
  ncdf->$VARNAME(ii)$ = SNOWFALL(:,:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
end if

delete(AGGRO_OUT_NAME)

end do

wallClockElapseTime(wcStrtWriteout, "Time to write output", 0)
wallClockElapseTime(wcStrt, "Total run time", 0)













end
