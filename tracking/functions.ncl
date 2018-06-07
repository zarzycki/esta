undef("getTrajectories")
function getTrajectories(tfile,ncol,maxLines)

; maxLines -- max # of lines allocated, but if set to negative integer, will try and
; determine internally. Only use hardcoded (e.g., 150) if you segfault or other OB errors

begin

;tfile="/glade/u/home/zarzycki/ncl/projects/slab/calcSST/trajectories.txt"
; read in trajectory data from ascii file
; get number of storms by counting number of lines with 'start'
nstorms_tot = 0
do ii = 0,120
  if (nstorms_tot .le. 1 .or. nstorms_tot .ge. 10000) then
    if (ii .ne. 0) then
      print("ii "+ii)
      system("sleep 2")
    end if
    nstorms_tot = stringtoint(systemfunc("grep 'start' "+tfile+" | wc -l"))
  else
    break
  end if
  print("GETTRAJECTORIES: Getting trajectories from: "+tfile)
  print("GETTRAJECTORIES: Total number of storms = "+nstorms_tot)
end do

tdata_tmp = asciiread(tfile,-1,"string")  ; read in trajectory data as strings
dims = dimsizes(tdata_tmp)    ; get dimensions
nrows = dims(0)               ; get number of rows
;ncol = 12                     ; assign number of columns
;maxLines = 200

delete([/dims/])

; ====DO NOT TOUCH THIS==============
tab_space="	 "
; ===================================

; =============== FIGURE OUT MAX NUMBERS OF STORM TIMES ========

if (maxLines .le. 0)
  print("GETTRAJECTORIES: No max lines defined, trying to figure out!")
  randDate = systemfunc("date +%s%N")
	ML_filetmpFile="supertempcolinhi"+randDate+".txt"
	system("grep -n 'start' "+tfile+" |cut -f1 -d: > "+ML_filetmpFile)
	ML_tmpStormRowsData = asciiread(ML_filetmpFile,-1,"integer")
	ML_tmpStormRowsDiff = (ML_tmpStormRowsData(1:dimsizes(ML_tmpStormRowsData)-1)-ML_tmpStormRowsData(0:dimsizes(ML_tmpStormRowsData)-2))-1
	ML_finalStormLine=systemfunc("grep 'start' "+tfile+" | tail -1")
	ML_finalStormSplit = str_split(ML_finalStormLine,tab_space)
	ML_lastStormSize = toint(ML_finalStormSplit(1))
	maxLines = max((/max(ML_tmpStormRowsDiff),ML_lastStormSize/))
	system("rm "+ML_filetmpFile)
end if

print("GETTRAJECTORIES: Allocating... "+maxLines+" max lines per storm")

; ===============================================================

; create new array that will hold all data except rows containing start information
; add column to contain storm number
tdata = new((/nstorms_tot,ncol,maxLines/),float)

; create new array to hold all header information for each storm
dur     = new(nstorms_tot,integer)
s_year  = new(nstorms_tot,integer)
s_month = new(nstorms_tot,integer)
s_day   = new(nstorms_tot,integer)
s_hour  = new(nstorms_tot,integer)

; loop over lines of trajectory data, convert string to float
nstorm = -1   ; assign storm number
n2 = 0       ; counter data without header lines

do n1 = 0,nrows-1  ; loop over all data with header lines
  first = str_get_cols(tdata_tmp(n1),0,0)   ; obtain first character of row
  splitStr = str_split(tdata_tmp(n1),tab_space)
  if (first .eq. "s") then              ; header line
    ;print("Header line")
    nstorm = nstorm + 1                 ; set storm number
    n2 = 0                              ; reset "per indiv traj" counter
    dur(nstorm)     = stringtoint(splitStr(1))
    s_year(nstorm)  = stringtoint(splitStr(2))
    s_month(nstorm) = stringtoint(splitStr(3))
    s_day(nstorm)   = stringtoint(splitStr(4))
    s_hour(nstorm)  = stringtoint(splitStr(5))
  else  ; if 's' in the header
    do z = 0,ncol-1
      ;print(nstorm+" "+n2+" "+z)
      ;print(" "+tdata_tmp(n1))
      tdata(nstorm,z,n2) = stringtofloat(splitStr(z))
    end do
    n2 = n2 + 1
  end if
  delete(splitStr)
end do

delete([/tdata_tmp,n1,n2,nstorms_tot,nstorm,first,nrows/])

return(tdata)

end



undef("bearing")
function bearing(
  lat1[1] : numeric,   ; lat1 in degrees
  lon1[1] : numeric,   ; lon1 in degrees
  lat2[1] : numeric,   ; lat2 in degrees
  lon2[1] : numeric    ; lon2 in degrees
  )
local pi, d2r, r2d, lat1r, lon1r, lat2r, lon2r, ang_tmp, lat1_tmp, lon1_tmp
; calculates the bearing (angle of movement over the sphere)
; between two points on a sphere
;
; output angle in degrees
begin

; constants
  pi = 4.*atan(1.)
  d2r = pi/180.
  r2d = (1/d2r)

; convert lat/lon points from degrees to radians
  lat1r = lat1*d2r
  lon1r = lon1*d2r
  lat2r = lat2*d2r
  lon2r = lon2*d2r

; compute angle of motion based on two lat/lon points
  ang_tmp = r2d*atan2(sin((lon2r-lon1r))*cos(lat2r), \
                  cos(lat1r)*sin(lat2r) - sin(lat1r)*cos(lat2r)*cos(lat2r-lat1r))

; convert angles from (-180 to +180) to (0 to 360)
  ang = mod(ang_tmp + 360.,360.)

  delete([/pi, d2r, r2d, lat1r, lon1r, lat2r, lon2r, ang_tmp/])

  return(ang)

end
