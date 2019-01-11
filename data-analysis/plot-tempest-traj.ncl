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

; User settings!
basin = "global"  ;Valid values: atlantic, wpac, tropics, (anything else for global)
output_type = "png"

ms_to_kts = 1.94384449

;thefile="tempest.LENS.5e-9_12.002"
thefile="trajectories.txt.LENS.5e-9_12.002"
ncols = 9 ; number of columns in non-header rows of trajectory file
maxLines = -1  ; set to -1 for auto-detection
traj_data = getTrajectories(thefile,ncols,maxLines)

; need to be modified to 2,3,4 if using latlon data
lon = traj_data(:,2,:)
lat = traj_data(:,3,:)
pres = traj_data(:,4,:)/100.

ndims = dimsizes(lon)
nstorms = ndims(0)
ntimes = ndims(1)

;--------------------------------------------------------
print("Beginning plot sequence")
wks  = gsn_open_wks(output_type,"latlon_traj_plot")
map = new(1,graphic)
  
res  = True
res@gsnDraw    = False
res@gsnFrame   = False
res@gsnMaximize = True
res@mpFillOn             = True              ; turn off gray continents
res@mpLandFillColor = "Tan"
res@mpOceanFillColor = "LightBlue1"
res@mpOutlineOn = True
res@mpInlandWaterFillColor = res@mpOceanFillColor
res@mpGreatCircleLinesOn = True
 
if (basin .eq. "atlantic")  
  res@mpMinLatF             = 10
  res@mpMaxLatF             = 70
  res@mpMinLonF             = 240
  res@mpMaxLonF             = 355
end if
if (basin .eq. "northeast")
  res@mpMinLatF             = 20
  res@mpMaxLatF             = 60
  res@mpMinLonF             = 270
  res@mpMaxLonF             = 330
end if
if (basin .eq. "nhemi")
  res@mpMinLatF             = 0
  res@mpMaxLatF             = 90
  res@mpMinLonF             = -150
  res@mpMaxLonF             = 120
end if
if (basin .eq. "wpac")  
  res@mpMinLatF             = 5
  res@mpMaxLatF             = 45
  res@mpMinLonF             = 115
  res@mpMaxLonF             = 180
end if
if (basin .eq. "tropics")
  res@mpMinLatF             = -65
  res@mpMaxLatF             = 65
  res@mpCenterLonF = 200.
end if
 
map = gsn_csm_map_ce(wks,res)

gsn_define_colormap(wks,"rainbow+white+gray")
  polyres               = True          ; poly marker mods desired 
  polyres@gsMarkerIndex = 1            ; choose circle as polymarker
  lineres = True





dum = new((/nstorms,ntimes/),graphic)

do i = 0,nstorms-1
do j = 0,ntimes-1

thisPres=pres(i,j)


linethicknesses = (/0.2,0.50,1.8,2.6,3.2,4.4/)


if (.not.ismissing(thisPres))

  if (thisPres .gt. 1005)
    polyres@gsMarkerColor = (/"blue"/)
    lineres@gsLineColor = (/"blue"/)
    lineres@gsLineThicknessF = linethicknesses(0)
  else if (thisPres .ge. 995 .and. thisPres .lt. 1005)
    polyres@gsMarkerColor = (/"green3"/)
    lineres@gsLineColor = (/"green3"/)
    lineres@gsLineThicknessF = linethicknesses(1)
  else if (thisPres .ge. 980 .and. thisPres .lt. 995)
    polyres@gsMarkerColor = (/"yellow2"/)
    lineres@gsLineColor = (/"yellow2"/)
    lineres@gsLineThicknessF = linethicknesses(2)
  else if (thisPres .ge. 965 .and. thisPres .lt. 980)
    polyres@gsMarkerColor = (/"orange"/)
    lineres@gsLineColor = (/"orange"/)
    lineres@gsLineThicknessF = linethicknesses(3)
  else if (thisPres .ge. 950 .and. thisPres .lt. 965)
    polyres@gsMarkerColor = (/"darkorange3"/)
    lineres@gsLineColor = (/"darkorange3"/)
    lineres@gsLineThicknessF  = linethicknesses(4)
  else
    polyres@gsMarkerColor = (/"red"/) 
    lineres@gsLineColor = (/"red"/)
    lineres@gsLineThicknessF  = linethicknesses(5)
  end if
  end if
  end if
  end if
  end if
  
 txres               = True                            
 txres@txFontHeightF = 0.01             ; Set the font height
  
  if (j .ne. 0) then
    dum(i,j) = gsn_add_polyline(wks,map,(/lon(i,j),lon(i,j-1)/),(/lat(i,j),lat(i,j-1)/),lineres)
  end if
  
end if

end do
end do





;******************
; Plot a legend   *
;******************

lgres                    = True

lgres@lgAutoManage       = False

lgres@vpWidthF           = 0.17       ; was 0.08        ; width of legend (NDC)
lgres@vpHeightF          = 0.13       ; was 0.08        ; height of legend (NDC)
; lgres@lgBottomMarginF    = 0.17     ; was 0.25

lgres@lgPerimFill        = 0                     ; Use solid fill (0) instead of the default hollow fill
lgres@lgPerimFillColor   = "Background"

;lgres@lgOrientation = "Horizontal"

; lgres@lgBoxMajorExtentF  = 0.4
lgres@lgBoxMinorExtentF  = 0.2	; controls how wide the box holding the legend items (lines and markers) can be in relation to legend

; lgres@lgBoxBackground    = "PaleTurquoise3"

lgres@lgMonoItemType        = False                 ; indicates that we wish to set the item types individually
lgres@lgMonoMarkerIndex     = False
lgres@lgMonoLineThickness   = False
lgres@lgMonoMarkerThickness = False
lgres@lgMonoMarkerSize      = False

;lgres@lgLabelFont             = 0
lgres@lgLabelFontHeightF      = 0.09
lgres@lgLabelFontAspectF      = 1.4
lgres@lgLabelConstantSpacingF = 0.0

legendType = "Markers"
markerSize = .01
markerIndex = 16

lgres@lgItemCount        = 6
;lgres@lgLineDashSegLenF  = 0.8
lgres@lgItemTypes        = (/"Markers","Markers","Markers","Markers","Markers","Markers"/)
lgres@lgMarkerIndexes    = (/      16,      16,      16,      16,      16,      16/)
lgres@lgMarkerThicknesses  = (/    markerSize,    markerSize,    markerSize,    markerSize,    markerSize,    markerSize/)
lgres@lgMarkerColors       = (/ "red",     "darkorange3",   "orange",  "yellow",   "green",    "blue"    /) ; colors for legend lines
;lgres@lgDashIndexes      = (/      0,          0,          0,         0,               0,      0      /) ; dash indexes
legend_labels = (/"<950mb","950-965mb","965-980mb",\
               "980-995mb","995-1005mb",">1005mb"/)

legend = gsn_create_legend(wks,lgres@lgItemCount,legend_labels,lgres)

amres = True
if (basin .eq. "atlantic") then
amres@amParallelPosF   = 0.43		    ; move legend to the right
amres@amOrthogonalPosF = -.3             ; move the legend down
else if (basin .eq. "tropics") then
amres@amParallelPosF   = 0.21		    ; move legend to the right .37 formerly
amres@amOrthogonalPosF = .27             ; move the legend down
else if (basin .eq. "nhemi") then
amres@amParallelPosF   = 0.30		    ; move legend to the right .37 formerly
amres@amOrthogonalPosF = .02             ; move the legend down
else
amres@amParallelPosF   = -0.30		    ; move legend to the right .37 formerly
amres@amOrthogonalPosF = .27             ; move the legend down
end if
end if
end if
annoid1 = gsn_add_annotation(map,legend,amres)   ; add legend to plot
  
print("Drawing plot")

resP = True
resP@amJust           = "TopLeft"

;  resP@gsnPanelFigureStrings= (/"a."/)

resP@gsnPanelFigureStringsFontHeightF = 0.027
gsn_panel(wks,map,(/1,1/),resP)             ; now draw as one plot

end
