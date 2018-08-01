;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/getTrajectories.ncl"

;========================================================
begin

wcStrt     = systemfunc("date")

; NL VARS
snowFileFull="/glade/u/home/zarzycki/scratch/LENS-snow/b.e11.BRCP85C5CNBDRD.f09_g16.002.cam.h2.PRECT_SNOW.2071010100Z-2080123118Z.nc"
SNOWVARNAME="PRECT_SNOW"
traj_filename="./traj/conntraj.LENS.2071.002.5e-9_12"
outDir="/glade/scratch/zarzycki/LES-snow/storm-files/"
outFileName="ind-storms.LENS.2071.002.5e-9_12.nc_TEST.nc"

doPrecip=False
doPSL=False

; How big of a window around track lat/lon do we want to take?
deltaLat=15.0
deltaLon=15.0

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

;-------------lookup table------------------------------

; Here we build a "lookup" table to find what times are on what files...
; get list of h3 files
filelist = systemfunc ("ls "+snowFileFull)
ftmp = addfile(filelist(0),"r")
nfiles = dimsizes(filelist)
nftimes = dimsizes(ftmp->time)
timeArr = new((/nfiles,nftimes/),typeof(ftmp->time))
do ii = 0,nfiles-1
;do ii = 0,11
  print("loading time for file #"+ii+" of "+(nfiles-1))
  ftmp = addfile(filelist(ii),"r")
  nftimes=dimsizes(ftmp->time)   ; needed since files at end of series may have less than rest
  timeArr(ii,0:nftimes-1) = ftmp->time
end do

; can get these from any arbitrary file
lat=ftmp->lat({minLatOut:maxLatOut})
lon=ftmp->lon({minLonOut:maxLonOut})
timeunits=ftmp->time@units
opt=0
opt@calendar=ftmp->time@calendar
delete(ftmp)

;-------------setup arrays------------------------------

wcStrtExtract     = systemfunc("date")

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

sttime=new(nstorms,"integer")
entime=new(nstorms,"integer")
durtime=new(nstorms,"integer")
sttime!0    = "storm"
entime!0    = "storm"
durtime!0   = "storm"

loaded_file=""   ; set loaded file to empty
do j = 0,nstorms-1
  print(" Storm #"+(j+1)+" of "+nstorms)
  durtime(j)=0
  do i = 0,ntimes-1
    if (.not.ismissing(traj_lon(j,i))) then

      print(i+"    "+ traj_yyyy(j,i)+" "+traj_mm(j,i)+" "+traj_dd(j,i)+" "+traj_hh(j,i) )
      thisTime=cd_inv_calendar(traj_yyyy(j,i),traj_mm(j,i),traj_dd(j,i),traj_hh(j,i),0.0,0.0,timeunits,opt)
      if (i .eq. 0) then
        sttime(j)=toint(traj_yyyy(j,i))*1000000+toint(traj_mm(j,i))*10000+toint(traj_dd(j,i))*100+toint(traj_hh(j,i))
      end if
      durtime(j) = durtime(j)+6
      entime(j)=toint(traj_yyyy(j,i))*1000000+toint(traj_mm(j,i))*10000+toint(traj_dd(j,i))*100+toint(traj_hh(j,i))
    
      ; look up time ix
      time1d      = ndtooned(timeArr)
      dsizes_time = dimsizes(timeArr)
      indices  = ind_resolve(ind(time1d.eq.thisTime),dsizes_time)
      fileix     = indices(0,0)
      trackindex = indices(0,1)
      needed_file=filelist(fileix)
      if (needed_file .ne. loaded_file) then
        ; need to load file
        print("Currently loading... "+needed_file)
        f1 = addfile(needed_file,"r")
        loaded_file=needed_file
      else
        print("existing loaded file... "+loaded_file)
      end if

      thislat=traj_lat(j,i)
      thislon=traj_lon(j,i)
      thispres=traj_pres(j,i)

      maxLat=thislat+deltaLat
      minLat=thislat-deltaLat
      maxLon=thislon+deltaLon
      minLon=thislon-deltaLon

      gcdist=gc_latlon(thislat,thislon,lat2D,lon2D,2,2)

      print(trackindex+"")
      SNOWFALL(j,i,:,:)=tofloat((f1->$SNOWVARNAME$(trackindex,{minLatOut:maxLatOut},{minLonOut:maxLonOut})))
      SNOWFALL(j,i,:,:)=where(gcdist.lt.deltaLat,SNOWFALL(j,i,:,:),SNOWFALL@_FillValue)
      SNOWFALL(j,i,:,:)=SNOWFALL(j,i,:,:)*6*ms_to_inhr
    end if
  end do
end do

CUM_SNOWFALL=dim_sum_n_Wrap(SNOWFALL,1)

wallClockElapseTime(wcStrtExtract, "Time to extract all snowfall events", 0)

; WRITE OUTPUT FILE

wcStrtWriteout     = systemfunc("date")

setfileoption("nc","Format","LargeFile")

outFile=outDir+"/"+outFileName

print("writing:   "+outFile)

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

; print 3-D vars
if (doTimeSeries) then
  ncdf->SNOWFALL = SNOWFALL(:,:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
end if

wallClockElapseTime(wcStrtWriteout, "Time to write output", 0)
wallClockElapseTime(wcStrt, "Total run time", 0)













end
