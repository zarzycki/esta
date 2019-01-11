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

;snowFileFull="/glade/u/home/zarzycki/scratch/LENS-snow/b.e11.B20TRC5CNBDRD.f09_g16.035.cam.h2.PRECT_SNOW.1990010100Z-2005123118Z.nc"
;SNOWVARNAME="PRECT_SNOW"
;traj_filename="/glade/u/home/zarzycki/snow-tracking/tracking/traj/tempest.LENS.1990.035"
;outDir="/glade/scratch/zarzycki/LES-snow/storm-files/"
;outFileName="test.nc"

; How big of a window around track lat/lon do we want to take?
;deltaLat=12.0
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

nlines=toint(systemfunc("wc -l < /glade/u/home/zarzycki/snow-tracking/tracking/reanalysis-track/filelist.SNOW.txt"))
all_files = asciiread("/glade/u/home/zarzycki/snow-tracking/tracking/reanalysis-track/filelist.SNOW.txt",nlines,"string")
fall = addfiles (all_files, "r")

;----------------------------------------------------------
; set up arrays

wcStrtExtract     = systemfunc("date")

lat=fall[0]->lat({minLatOut:maxLatOut})
lon=fall[0]->lon({minLonOut:maxLonOut})
nlat=dimsizes(lat)
nlon=dimsizes(lon)
lat2D=conform_dims((/nlat,nlon/),lat,0) 
lon2D=conform_dims((/nlat,nlon/),lon,1) 

sttime=new(nstorms,"integer")
entime=new(nstorms,"integer")
durtime=new(nstorms,"integer")
sttime!0    = "storm"
entime!0    = "storm"
durtime!0   = "storm"

;----------------------------------------------------------

PREC_FILE=addfile(snowFileFull,"r")

filetime=fall[:]->time
maxFileTime=max(filetime)
minFileTime=min(filetime)
FULLPREC_DOUBLE=fall[:]->$SNOWVARNAME$(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
FULLPREC=tofloat(FULLPREC_DOUBLE)
;FULLPREC=ddd
copy_VarMeta(FULLPREC_DOUBLE,FULLPREC)
delete(FULLPREC_DOUBLE)
STORMPREC=FULLPREC
STORMHITS=FULLPREC
STORMHITS=0.0
printVarSummary(FULLPREC)

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
  
  if (.not. isdefined("units"))
    units=fall[0]->time@units
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
    STORMPREC(trackindex,:,:)=where(gcdist.ge.deltaLat,FULLPREC(trackindex,:,:),0.0)
    STORMHITS(trackindex,:,:)=where(gcdist.lt.deltaLat,STORMHITS(trackindex,:,:)+1,STORMHITS(trackindex,:,:))
  end if   ; check to see if tracktime on file
  
end if

end do
end do

DOUBLESNOW=FULLPREC
DOUBLESNOW=where(STORMHITS.gt.1.5,FULLPREC,0.0)

NONSTORMSUM=dim_sum_n_Wrap(STORMPREC,0)
TOTALSUM=dim_sum_n_Wrap(FULLPREC,0)
DOUBLESUM=dim_sum_n_Wrap(DOUBLESNOW,0)
STORMSUM=NONSTORMSUM
STORMSUM=TOTALSUM-NONSTORMSUM
TOTALSUM=where(TOTALSUM.lt.1e-8,TOTALSUM@_FillValue,TOTALSUM)

STORMFRAC = TOTALSUM
STORMFRAC = STORMSUM / TOTALSUM

DOUBLEFRAC = TOTALSUM
DOUBLEFRAC = DOUBLESUM / TOTALSUM

wallClockElapseTime(wcStrtExtract, "Time to extract all snowfall events", 0)


if (0 .eq. 1) then
  ; PLOT FRAC
  wks = gsn_open_wks("x11","snowfraction")               ; open workstation

  res                     = True         ; plot modifications desired
  res@gsnDraw=False
  res@gsnFrame=False

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

  res@mpOceanFillColor = "White"
  res@cnFillDrawOrder             = "PreDraw"
  res@cnLineDrawOrder             = "PreDraw"
  res@mpFillAreaSpecifiers        = (/"Land"/)
  res@mpSpecifiedFillColors       = (/"transparent"/)
  res@mpOutlineDrawOrder = "Draw"

  res@cnFillPalette        = "WhiteBlueGreenYellowRed"              ; choose colormap
  res@cnSpanFillPalette=True

  ;res@cnFillMode          = "RasterFill"
  res@gsnAddCyclic = False

  res@mpMinLatF             = 38.
  res@mpMaxLatF             = 48.
  res@mpMinLonF             = 279.
  res@mpMaxLonF             = 293.
  res@mpCenterLonF = 286.

  res@mpOutlineOn                 = True

  res@mpUSStateLineThicknessF     = 2.5  ; 2-1/2 times as thick.
  res@mpGeophysicalLineThicknessF = 2.5

  res@mpOutlineBoundarySets       = "GeophysicalAndUSStates"
  res@mpOutlineSpecifiers         = (/"Land","Conterminous US:states"/)

  res@gsnStringFontHeightF=0.027

  resavg=res

  resavg@tiMainString=""
  resavg@gsnRightString=""
  resavg@gsnLeftString="Fraction of snowfall"
  map = gsn_csm_contour_map_ce(wks,STORMFRAC,resavg)
  print(""+min(STORMFRAC))
  print(""+max(STORMFRAC))

  draw(map)
  frame(wks)
end if

; WRITE OUTPUT FILE

wcStrtWriteout     = systemfunc("date")

setfileoption("nc","Format","LargeFile")

outFile=outDir+"/"+outFileName

print("Writing... "+outFile)

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
;ncdf->STORMPREC = STORMPREC(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
ncdf->STORMHITS = STORMHITS(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
;ncdf->FULLPREC = FULLPREC(:,{minLatOut:maxLatOut},{minLonOut:maxLonOut})
ncdf->STORMFRAC = STORMFRAC({minLatOut:maxLatOut},{minLonOut:maxLonOut})
ncdf->DOUBLEFRAC = DOUBLEFRAC({minLatOut:maxLatOut},{minLonOut:maxLonOut})

wallClockElapseTime(wcStrtWriteout, "Time to write output", 0)
wallClockElapseTime(wcStrt, "Total run time", 0)

end
