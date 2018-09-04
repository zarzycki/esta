;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/surfaceWindCorrFactor.ncl"
load "/glade/u/home/zarzycki/ncl/functions/getTrajectories.ncl"
;========================================================
begin

wcStrt     = systemfunc("date")

; Check for vars from command line
if (any(isvar((/"prectFileFull","PRECTVARNAME"/)))) then
  if (all(isvar((/"prectFileFull","PRECTVARNAME"/)))) then
    doPrecip=True
    print("doPrecip set to true")
  else
    print("only one of prectFileFull or PRECTVARNAME is defined!")
    exit
  end if
end if
if (any(isvar((/"pslFileFull","PSLVARNAME"/)))) then
  if (all(isvar((/"pslFileFull","PSLVARNAME"/)))) then
    doPSL=True
    print("doPSL set to true")
  else
    print("only one of pslFileFull or PSLVARNAME is defined!")
    exit
  end if
end if

; How big of a window around track lat/lon do we want to take?
;deltaLat=10.0
print("read in deltaLat: "+deltaLat+" degrees.")
deltaLon=deltaLat

doTimeSeries=True
doRegionalOutput=True

; If yes, truncation output
if (doRegionalOutput) then
  minLatOut=22.
  maxLatOut=55.
  minLonOut=265.
  maxLonOut=307.
else
  minLatOut=-90.
  maxLatOut=90.
  minLonOut=0.
  maxLonOut=360.
end if

;----------------Constants-----------------------------------

ms_to_kts = 1.94384449
pi = 3.14159265358979
omega = 7.2921E-5
ms_to_inhr=141732.

;-------------Get trajectories------------------------------

thefile=traj_filename
ncols = 9                     ; assign number of columns
maxLines = -1
traj_data = getTrajectories(thefile,ncols,maxLines)

traj_lon = traj_data(:,2,:)
traj_lat = traj_data(:,3,:)
traj_pres = traj_data(:,4,:)/100.
traj_yyyy = traj_data(:,5,:)
traj_mm = traj_data(:,6,:)
traj_dd = traj_data(:,7,:) 
traj_hh = traj_data(:,8,:)

ndims = dimsizes(traj_lon)
nstorms = ndims(0)
storm_index=ispan(0,nstorms-1,1)
ntimes = ndims(1)

;----------------------------------------------------------
; set up arrays

wcStrtExtract     = systemfunc("date")

fileWLATLON=addfile(snowFileFull,"r")
lat=fileWLATLON->lat({minLatOut:maxLatOut})
lon=fileWLATLON->lon({minLonOut:maxLonOut})
nlat=dimsizes(lat)
nlon=dimsizes(lon)
lat2D=conform_dims((/nlat,nlon/),lat,0) 
lon2D=conform_dims((/nlat,nlon/),lon,1) 
varArray=new((/nstorms,ntimes,nlat,nlon/),"double")

varArray!0    = "storm"
varArray!1    = "time"
varArray!2    = "lat"
varArray!3    = "lon"
varArray&lat = lat
varArray&lon = lon

SNOWFALL=varArray
if (doPrecip) then
  PRECT=varArray
end if
if (doPSL) then
  PSL=varArray
end if

sttime=new(nstorms,"integer")
entime=new(nstorms,"integer")
durtime=new(nstorms,"integer")
sttime!0    = "storm"
entime!0    = "storm"
durtime!0   = "storm"

;----------------------------------------------------------

PREC_FILE=addfile(snowFileFull,"r")
if (doPrecip) then
  PRECT_FILE=addfile(prectFileFull,"r")
end if
if (doPSL) then
  PSL_FILE=addfile(pslFileFull,"r")
end if

filetime=PREC_FILE->time
maxFileTime=max(filetime)
minFileTime=min(filetime)

opt=0   ; cd_inv_calendar
; check to see if calendar attribute is attached, if so, use for cd_inv_calendar
if ( isatt(filetime,"calendar") ) then
  print("Using calendar: "+filetime@calendar)
  opt@calendar=filetime@calendar
end if

do j = 0,nstorms-1
print(" Storm #"+(j+1)+" of "+nstorms)
durtime(j)=0
do i = 0,ntimes-1
if (.not.ismissing(traj_lon(j,i))) then
  thisYear=traj_yyyy(j,i)
  ;units="days since 2006-01-01 00:00:00"
  ;units="days since 1920-01-01 00:00:00"
  if (.not. isdefined("units"))
    units=PREC_FILE->time@units
  end if

  yyyy=toint(traj_yyyy(j,i))
  mm=toint(traj_mm(j,i))
  dd=toint(traj_dd(j,i))
  hh=toint(traj_hh(j,i))
  sssss=hh*3600

  print(yyyy+" "+mm+" "+dd+" "+hh+" "+sssss)

  if (i .eq. 0) then
    sttime(j)=yyyy*1000000+mm*10000+dd*100+hh
  end if
  durtime(j) = durtime(j)+6
  entime(j)=yyyy*1000000+mm*10000+dd*100+hh

  tracktime = cd_inv_calendar(yyyy,mm,dd,hh,0,0,units,opt)
  trackindex = ind(filetime .eq. tracktime)

  thislat=traj_lat(j,i)
  thislon=traj_lon(j,i)
  thispres=traj_pres(j,i)

  maxLat=thislat+deltaLat
  minLat=thislat-deltaLat
  maxLon=thislon+deltaLon
  minLon=thislon-deltaLon

  gcdist=gc_latlon(thislat,thislon,lat2D,lon2D,2,2)

  if ( tracktime .le. maxFileTime .and. tracktime .ge. minFileTime ) then
    print(trackindex+"")
    SNOWFALL(j,i,:,:)=tofloat((PREC_FILE->$SNOWVARNAME$(trackindex,{minLatOut:maxLatOut},{minLonOut:maxLonOut})))
    SNOWFALL(j,i,:,:)=where(gcdist.lt.deltaLat,SNOWFALL(j,i,:,:),SNOWFALL@_FillValue)
    ;SNOWFALL(j,i,:,:)=where(lat2D.le.maxLat .and. lat2D.ge.minLat,SNOWFALL(j,i,:,:),SNOWFALL@_FillValue)
    ;SNOWFALL(j,i,:,:)=where(lon2D.le.maxLon .and. lon2D.ge.minLon,SNOWFALL(j,i,:,:),SNOWFALL@_FillValue)
    SNOWFALL(j,i,:,:)=SNOWFALL(j,i,:,:)*6*ms_to_inhr
    if (doPrecip) then
      PRECT(j,i,:,:)=tofloat((PRECT_FILE->PRECT(trackindex,{minLatOut:maxLatOut},{minLonOut:maxLonOut})))
      PRECT(j,i,:,:)=where(gcdist.lt.deltaLat,PRECT(j,i,:,:),PRECT@_FillValue)
      ;PRECT(j,i,:,:)=where(lat2D.le.maxLat .and. lat2D.ge.minLat,PRECT(j,i,:,:),PRECT@_FillValue)
      ;PRECT(j,i,:,:)=where(lon2D.le.maxLon .and. lon2D.ge.minLon,PRECT(j,i,:,:),PRECT@_FillValue)
      PRECT(j,i,:,:)=PRECT(j,i,:,:)*6*ms_to_inhr
    end if
    if (doPSL) then
      PSL(j,i,:,:)=tofloat((PSL_FILE->PSL(trackindex,{minLatOut:maxLatOut},{minLonOut:maxLonOut})))
      ;PSL(j,i,:,:)=where(lat2D.le.maxLat .and. lat2D.ge.minLat,PSL(j,i,:,:),PSL@_FillValue)
      ;PSL(j,i,:,:)=where(lon2D.le.maxLon .and. lon2D.ge.minLon,PSL(j,i,:,:),PSL@_FillValue)
      PSL(j,i,:,:)=where(gcdist.lt.deltaLat,PSL(j,i,:,:),PSL@_FillValue)
    end if
  end if   ; check to see if tracktime on file
  
end if

end do
end do

CUM_SNOWFALL=dim_sum_n_Wrap(SNOWFALL,1)
if (doPrecip) then
  CUM_PRECT=dim_sum_n_Wrap(PRECT,1)
end if
if (doPSL) then
  MIN_PSL=dim_min_n_Wrap(PSL,1)
end if

wallClockElapseTime(wcStrtExtract, "Time to extract all snowfall events", 0)

; WRITE OUTPUT FILE

wcStrtWriteout     = systemfunc("date")

setfileoption("nc","Format","LargeFile")

outFile=outDir+"/"+outFileName

system("mkdir -p "+outDir)   ; remove any pre-existing file
system("/bin/rm -f "+outFile)   ; remove any pre-existing file
ncdf = addfile(outFile,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

; print 1-D vars
ncdf->sttime=sttime
ncdf->entime=entime
ncdf->durtime=durtime

; print 2-D vars
ncdf->CUM_SNOWFALL = CUM_SNOWFALL(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
if (doPrecip) then
  ncdf->CUM_PRECT = CUM_PRECT(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
end if
if (doPSL) then
  ncdf->MIN_PSL = MIN_PSL(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
end if

; print 3-D vars
if (doTimeSeries) then
  ncdf->SNOWFALL = SNOWFALL(:,:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
  if (doPrecip) then
    ncdf->PRECT = PRECT(:,:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
  end if
  if (doPSL) then
    ncdf->PSL = PSL(:,:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
  end if
end if

wallClockElapseTime(wcStrtWriteout, "Time to write output", 0)
wallClockElapseTime(wcStrt, "Total run time", 0)

end
