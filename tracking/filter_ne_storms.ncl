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
fcntoload=ESTAPATH+"/functions/functions.ncl"
if (fileexists(fcntoload)) then
  loadscript(fcntoload)
else
  print("function "+fcntoload+" doesn't exist")
  exit
end if

;========================================================
begin

;infile="trajectories.txt.MP15A-120A-US-EXP213"
;outfile="traj.filt"
;threshold = 20.0

system("rm "+outfile)

searchLat = 40.
searchLon1 = 288.
searchLon2 = 285.
threshold = 8.0
; usually 10 and 4
firstmm=10
lastmm=4

; need to be modified to 2,3,4,6 if using LENS/latlon data
thefile=infile
ncols = 8                     ; assign number of columns
maxLines = -1
headerDelimStr = "start"
delimStr="tab"
isUnstruc=False
isHeader=False
traj_data = getTrajectoriesGeneric(infile,ncols,maxLines,headerDelimStr,delimStr,isUnstruc,isHeader)

lon = traj_data(:,2,:)
lat = traj_data(:,3,:)
pres = traj_data(:,4,:)/100.
year = traj_data(:,5,:)
month = traj_data(:,6,:)
ndims = dimsizes(lon)
nstorms = ndims(0)
ntimesall = ndims(1)

stcounter=1
encounter=-1

;sed -n '8,12p' yourfile

do ii = 0,nstorms-1
  print("Doing "+ii+" of: "+(nstorms-1))
  onedpres = pres(ii,:)
  onedlon = lon(ii,:)
  onedlat = lat(ii,:)
  onedmonth = month(ii,:)
  onedyear = year(ii,:)
  ntimes = num(.not.(ismissing(onedlon)))
  onedbearing = onedlon
  onedbearing(0) = onedbearing@_FillValue
  do k = 1,ntimes-1
    onedbearing(k) = bearing(onedlat(k-1),onedlon(k-1),onedlat(k),onedlon(k))
  end do

  if (ii .gt. 0)
    stcounter=encounter+1
  end if
  encounter = stcounter + ntimes

  isStorm=new(dimsizes(onedlon),"integer")
  isStorm=1
  gc1=gc_latlon(onedlat,onedlon,conform(onedlat,searchLat,-1),conform(onedlat,searchLon1,-1),2,2)
  gc2=gc_latlon(onedlat,onedlon,conform(onedlat,searchLat,-1),conform(onedlat,searchLon2,-1),2,2)
  isStorm=where(ismissing(onedbearing),isStorm@_FillValue,isStorm)
  isStorm=where(gc1 .le. threshold .OR. gc2 .le. threshold,isStorm,isStorm@_FillValue)
  ;isStorm=where(onedbearing .ge. 0 .AND. onedbearing .le. 90 ,isStorm,isStorm@_FillValue)

  ;if ((onedmonth(0) .ge. firstmm .or. onedmonth(0) .le. lastmm) .and. onedyear(0) .lt. 2000)
  if ((onedmonth(0) .ge. firstmm .or. onedmonth(0) .le. lastmm))
    inSeason=True
  else
    inSeason=False
  end if

  if(num(.not.ismissing(isStorm)) .gt. 0 .and. inSeason .eq. True )
    print("*************** Keep trajectory!")
    system("sed -n '"+toint(stcounter)+","+toint(encounter)+"p' "+infile+" >> "+outfile)
  else
    print("Delete trajectory!")
  end if

  delete(gc1)
  delete(gc2)
  delete(onedpres)
  delete(onedlat)
  delete(onedlon)
  delete(onedmonth)
  delete(onedbearing)

end do

exit

end
