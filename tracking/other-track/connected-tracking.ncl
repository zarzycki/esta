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

;ensmember="021"
;timestring="1990010100Z-2005123118Z"
;snowFileFull="/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/PRECT/b.e11.B20TRC5CNBDRD.f09_g16."+ensmember+".cam.h2.PRECT.1990010100Z-2005123118Z.nc"

;snowFileFull="/glade/u/home/zarzycki/scratch/h1files/JRA/2000/JRA.h1.2000.PRECT.nc"
;outfile="conntraj.JRA.2000"

adjustLon = True

minLat=37.0
maxLat=48.0
minLon=278.75
maxLon=294.0

thresh_prec=tofloat(thresh_prec_str)
print("using a threshold of: "+thresh_prec+" m/s")

;thresh_prec=7.e-9
;thresh_prec=10.e-9

PREC_FILE=addfile(snowFileFull,"r")
lat=tofloat(PREC_FILE->lat({minLat:maxLat}))
lon=tofloat(PREC_FILE->lon({minLon:maxLon}))
time=PREC_FILE->time

deltaTime = time(1)-time(0)

; calculate weighting
pi = 3.14159
gw = cos(lat*pi/180.)

maskedTime=cd_calendar(time,-3)

PRECT=tofloat(PREC_FILE->PRECT(:,{minLat:maxLat},{minLon:maxLon}))
dims=dimsizes(PRECT)
ntimes=dims(0)
nlat=dims(1)
nlon=dims(2)

print(ntimes)

maskFileName="/glade/u/home/zarzycki/snow-tracking/calc_RSI/NEmask.nc"
maskFile=addfile(maskFileName,"r")
MASKNAT=maskFile->MASK
MASKNAT_lat=MASKNAT&lat
if (adjustLon) then
  MASKNAT_lon=MASKNAT&lon+360.
  MASKINT=linint2_Wrap(MASKNAT&lon+360.,MASKNAT&lat,MASKNAT,False,lon,lat, 0)
else
  MASKNAT_lon=MASKNAT&lon
  MASKINT=linint2_Wrap(MASKNAT&lon,MASKNAT&lat,MASKNAT,False,lon,lat, 0)
end if

PRECTINT=linint2_Wrap(lon,lat,PRECT,False,MASKNAT_lon,MASKNAT_lat, 0)

gauswt=conform_dims((/ntimes,nlat,nlon/),gw,1)
MASK=conform_dims((/ntimes,nlat,nlon/),MASKINT,(/1,2/))

printVarSummary(MASK)
printVarSummary(gauswt)
printVarSummary(PRECT)

printVarSummary(PRECTINT)
printVarSummary(MASKNAT)
MASKNAT_3D=conform_dims((/ntimes,dimsizes(MASKNAT_lat),dimsizes(MASKNAT_lon)/),MASKNAT,(/1,2/))
MASKNAT_3D=where(MASKNAT_3D.le.0.5,MASKNAT_3D@_FillValue,MASKNAT_3D)

MASKPREC=PRECTINT
MASKPREC=MASKNAT_3D*PRECTINT
MASKPREC=MASKPREC*3.402e6   ; inches/day

; trigger 1 if "precip event," 0 if not.
MASKPREC_TF=MASKPREC
MASKPREC_TF=0.0

MEANPREC=dim_sum_n_Wrap(PRECT,(/1,2/))
MEANPREC=dim_sum_n_Wrap(MASK*gauswt*PRECT,(/1,2/))/dim_sum_n_Wrap(MASK*gauswt,(/1,2/))

prec_logic=new(ntimes,"string")
prec_logic="FALSE"

do ii = 0,ntimes-1
  if (MEANPREC(ii) .gt. thresh_prec) then
    prec_logic(ii)="TRUE"
    MASKPREC_TF(ii,:,:)=1.0
  end if
  ;print(MEANPREC(ii)+" "+maskedTime(ii)+" "+prec_logic(ii))
end do 

; ===== build connectivity
storm_index = -1
storm_initialized=False
storm_start = new(90000,"double")
storm_end = new(90000,"double")

do jj = 0,ntimes-1
  leadString=maskedTime(jj)+" "+prec_logic(jj)+" "

  if (jj .eq. (ntimes-1))
    if (prec_logic(jj) .eq. "TRUE" .and. storm_initialized .eq. True)
      print(leadString+"storm "+storm_index+" terminated at final time: "+maskedTime(jj))
      storm_initialized = False
      storm_end(storm_index) = time(jj)
    ; we need this logic because otherwise code will try to search broken code trajectory below and go OOB.
    else if (prec_logic(jj) .eq. "FALSE" .and. storm_initialized .eq. True)
      print(leadString+"storm "+storm_index+" terminated at second-to-last time: "+maskedTime(jj-1))
      storm_initialized = False
      storm_end(storm_index) = time(jj-1)
    else
      print("Done tracking...")
    end if
    end if
    break  ; break out of loop
  end if

  if (prec_logic(jj) .eq. "TRUE" .and. storm_initialized .eq. False)
    ; check to make sure next timestep is true
    if (prec_logic(jj+1) .eq. "TRUE") then
      storm_initialized = True
      storm_index = storm_index+1
      storm_start(storm_index) = time(jj)
      print(leadString+"storm initialized as storm: "+storm_index)
    else
      print(leadString+"false alarm")
    end if
  else if (prec_logic(jj) .eq. "TRUE" .and. storm_initialized .eq. True)
    print(leadString+"continuing storm: "+storm_index)
  else if (prec_logic(jj) .eq. "FALSE" .and. storm_initialized .eq. False)
    print(leadString+"no storm, doing nothing.")
  else if (prec_logic(jj) .eq. "FALSE" .and. storm_initialized .eq. True)
    ; check next time for broken traj
    if (prec_logic(jj+1) .eq. "TRUE") then
      print(leadString+"broken trajectory, continuing storm: "+storm_index)
    else
      print(leadString+"storm "+storm_index+" terminated at previous time: "+maskedTime(jj-1))
      storm_initialized = False
      storm_end(storm_index) = time(jj-1)
    end if
  else
    print("Bad logic... exiting...")
  end if
  end if
  end if
  end if
end do

move_frac=deltaTime  ; 0.25 for days, 6.0 for hours

system("rm -f "+outfile)
quote = inttochar(34)
do kk = 0,storm_index
  ;print("Storm: "+kk+"  start: "+storm_start(kk)+"  end: "+storm_end(kk))
  storm_length=toint((storm_end(kk)-storm_start(kk))/move_frac)+1
  utc_date = cd_calendar(storm_start(kk), 0)
  if (toint(utc_date(:,1)) .ge. 10 .or. toint(utc_date(:,1)) .le. 4 .and. .not.ismissing(storm_end(kk))) then
    outString="start   "+storm_length+"       "+utc_date(:,0)+"    "+utc_date(:,1)+"       "+utc_date(:,2)+"       "+utc_date(:,3)
    print(""+outString)
    system("echo "+quote+outString+quote+" >> "+outfile)
    do ii = 0,storm_length-1
      thisTime=storm_start(kk)
      thisTime=thisTime+move_frac*ii
      utc_date = cd_calendar(thisTime, 0)  
      outString="        999     999     287.000000      43.000000      -999.99    "+utc_date(:,0)+"    "+utc_date(:,1)+"       "+utc_date(:,2)+"       "+utc_date(:,3)
      print(""+outString)
      system("echo "+quote+outString+quote+" >> "+outfile)
    end do
  end if
end do

exit

;-------------Get trajectories------------------------------

if (0 .eq. 1) then
  thefile=traj_dir+"/"+traj_filename
  ncols = 9                     ; assign number of columns
  maxLines = -1
  traj_data = getTrajectories(thefile,ncols,maxLines)

  traj_lon = traj_data(:,2,:)
  traj_lat = traj_data(:,3,:)
  traj_pres = traj_data(:,4,:)/100.
  traj_yyyy = toint(traj_data(:,5,:))
  traj_mm = toint(traj_data(:,6,:))
  traj_dd = toint(traj_data(:,7,:))
  traj_hh = toint(traj_data(:,8,:))


  trackedTime2D=traj_yyyy*1000000+traj_mm*10000+traj_dd*100+traj_hh
  trackedTime=toint(ndtooned(trackedTime2D))


  ;do jj = 0,ntimes-1
  ;  if (prec_logic(jj) .eq. "TRUE")
  ;    match=any(trackedTime.eq.maskedTime(jj))
  ;    if(match)
  ;      thestring="YES"
  ;    else
  ;      thestring="NO"
  ;    end if
  ;    print(maskedTime(jj)+" "+thestring)
  ;  end if
  ;end do

end if

system("/bin/rm -f /glade/scratch/zarzycki/simple.nc")   ; remove any pre-existing file
ncdf = addfile("/glade/scratch/zarzycki/simple.nc" ,"c")  ; open output netCDF file
fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Simple Approach to netCDF Creation"
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True) 

ncdf->MASKINT  = MASKINT   
ncdf->MASKPREC = MASKPREC
ncdf->MASKPREC_TF = MASKPREC_TF

end
