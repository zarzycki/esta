f = addfile("simple.nc","r")
f2 = addfile("../gen_pop_files/population.nc","r")

maskdata = f->maskdata
POP=f2->POP
POP=POP/1000.

printVarSummary(POP)
printVarSummary(maskdata)
print(min(maskdata)+" "+max(maskdata))

maskdata = where(maskdata.eq.0,maskdata@_FillValue,maskdata)

wks  = gsn_open_wks("x11","shapefiles") ; send graphics to PNG file

res                     = True

res@gsnDraw             = False       ; don't draw yet
res@gsnFrame            = False       ; don't advance frame yet
res@gsnMaximize         = True        ; maximize plot in frame

res@mpProjection        = "LambertConformal"    ; choose projection
res@mpLambertParallel1F = 33          ; first parallel
res@mpLambertParallel2F = 45          ; second parallel
res@mpLambertMeridianF  = -98         ; meridian

res@mpLimitMode         = "Corners"   ; corner method of zoom
res@mpLeftCornerLatF    = 22          ; left corner
res@mpLeftCornerLonF    = -125        ; left corner
res@mpRightCornerLatF   = 50          ; right corner
res@mpRightCornerLonF   = -64         ; right corner

res@pmTickMarkDisplayMode = "Always"  ; turn on tickmarks

res@tiMainString          = "Percentage unemployment, by state"

res@mpOutlineBoundarySets     = "geophysicalandusstates"; turn on states
res@mpDataBaseVersion         = "mediumres"             ; select database
res@mpDataSetName             = "Earth..2"

res@cnFillOn         = True               ; color Fill 
res@cnFillMode       = "RasterFill"       ; Raster Mode
res@cnLinesOn        =  False             ; Turn off contour lines

res@cnFillOpacityF=0.4

plot = gsn_csm_contour_map(wks,maskdata,res)   ; Create map, but don't draw it yet.

;

res@lbLabelBarOn = False
res@cnFillOpacityF=0.4

res@cnFillOn            = True         ; Turn on contour fill

res@cnLinesOn           = False        ; Turn off contour lines
res@cnLineLabelsOn      = False        ; Turn off contour line labels
res@cnInfoLabelOn = False

res@cnFillMode          = "RasterFill"
res@cnLevelSelectionMode = "ExplicitLevels"                     ; set explicit contour levels
cmap = read_colormap_file("GMT_gray")
res@cnFillPalette = cmap(::-1,:)     ; reverse color map
delete(cmap)

res@cnLevels    = (/ 10,15,25,50,100/) ; set levels

map = gsn_csm_contour(wks,POP,res)

overlay(plot,map)


draw(plot)
frame(wks)