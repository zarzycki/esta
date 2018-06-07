;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;========================================================
begin

f1=addfile("~/scratch/LENS-snow/b.e11.B20TRC5CNBDRD.f09_g16.022.cam.h2.PRECT_SNOW.1990010100Z-2005123118Z.nc","r")

lat=f1->lat
lon=f1->lon
nlat=dimsizes(lat)
nlon=dimsizes(lon)
lat2D=conform_dims((/nlat,nlon/),lat,0) 
lon2D=conform_dims((/nlat,nlon/),lon,1)

thislat=40.
thislon=-70.

gcdist=gc_latlon(thislat,thislon,lat2D,lon2D,2,2)

PRECT=f1->PRECT_SNOW(0,:,:)

PRECT=where(gcdist.lt.15.0,PRECT,PRECT@_FillValue)

wks = gsn_open_wks("x11","proj")          ; send graphics to PNG file

res                   = True
res@mpFillOn          = False
res@cnFillOn          = True              ; color plot desired
res@cnFillPalette     = "gui_default"     ; set color map
res@cnLineLabelsOn    = False             ; turn off contour lines

contour = gsn_csm_contour_map(wks,PRECT,res)  ; create the plot

end
