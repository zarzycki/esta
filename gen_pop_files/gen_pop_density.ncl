;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/calc_area.ncl"
;========================================================
begin

; SEDAC population (year 2000)
; http://neo.sci.gsfc.nasa.gov/view.php?datasetId=SEDAC_POP
; Got 0.1deg CSV from above website
;ncols         8688     ;362
;nrows         1320      ;55
;xllcorner     -181
;yllcorner     18
;cellsize      0.0416666666667  ;
;NODATA_value  -9999
;pop=asciiread("usads00ag.asc",(/1320,8688/),"float")

doplot=False

filename = "SEDAC_POP_2000-01-01_rgb_3600x1800.CSV"

;---Read the values in as 1D, since we don't know rows and columns yet.
values_1d = asciiread(filename,-1,"string")
printVarSummary(values_1d)
ncols     = dimsizes(str_split(values_1d(0),","))
nrows     = dimsizes(values_1d)

;thevals=str_split(values_1d,",")

;---Reshape as 2D array, if needed, and convert to integer.
;values_2d = tofloat(onedtond(values_1d,(/nrows,ncols/)))

print(ncols)
print(nrows)
values_2d = new((/nrows,ncols/),"float",99999.0) ;99999.0 is missing value from data, do not change!
do i = 0,nrows-1
  values_2d((nrows-1)-i,:) = tofloat(str_split(values_1d(i),","))
end do


csvLon=fspan(-179.95,179.95,ncols)
csvLat=fspan(-89.95,89.95,nrows)
csvLon@units="degrees_east"
csvLat@units="degrees_north"
values_2d!0="lat"
values_2d!1="lon"
values_2d&lat=csvLat
values_2d&lon=csvLon

area=calc_area(csvLat,csvLon)
area@_FillValue=-9.e7
print(sum(area))
copy_VarCoords(values_2d,area)

totPop=values_2d*area
copy_VarCoords(values_2d,totPop)
print(sum(totPop))

;=====================================================================
; Assume variables T, PS exist and that they have coordinate variables 
;      T  is four  dimensional:  T(time,lev,lat,lon)
;      PS is three dimensional: PS(time,lat,lon)
;
; Assume that ORO exists. ORO is two dimensional but it does not have any 
; meta data associated with it.  These must be explicitly assigned. Assume the
; lat/lon coordinate arrays can be taken from variable T. 
; Note that the variable name on the output file does not 
; have to be same as variable name in the code.  
;=====================================================================

print ("Generating netCDF file...")

system("/bin/rm -f population.nc")   ; remove any pre-existing file
ncdf = addfile("population.nc" ,"c")  ; open output netCDF file

; define file attributes
fAtt               = True            ; assign file attributes
fAtt@title         = "SEDAC population (year 2000)"
fAtt@source_file   =  filename
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

; write output vars
ncdf->DENS  = values_2d       
ncdf->AREA = area
ncdf->POP = totPop  

print ("... done")

if (doplot) then
  ;******************************************************
  ; create plot
  ;******************************************************
  wks = gsn_open_wks("x11","homme")               ; open workstation
  gsn_define_colormap(wks,"MPL_BuPu")   ; choose colormap

  res                     = True         ; plot modifications desired
  res@gsnSpreadColors     = True         ; Use full colormap, but start

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

  res@cnFillMode          = "RasterFill"
  res@gsnAddCyclic = False

  res@mpMinLatF             = 30.
  res@mpMaxLatF             = 50.
  res@mpMinLonF             = 265.
  res@mpMaxLonF             = 295.

  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels 
  res@cnLevels = ((/ 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000/)) ; set levels 

  res@mpOutlineOn                 = True

  res@mpUSStateLineThicknessF     = 1  ; 2-1/2 times as thick.
  res@mpGeophysicalLineThicknessF = 1

  res@mpOutlineBoundarySets       = "GeophysicalAndUSStates"
  res@mpOutlineSpecifiers         = (/"Land","Conterminous US:states"/)
  ;res@mpMaskOutlineSpecifiers     = (/"water"/)

  map = gsn_csm_contour_map_ce(wks,values_2d(:,:),res)

end if

end
