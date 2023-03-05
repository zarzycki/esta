load "./shapefile_utils.ncl"


masktype="NCA"

lonstart=-180.
gridsize=1.00

latS = -90.
latN = 90.
lonW = lonstart
lonE = lonstart+360.

dlat =  gridsize
dlon =  gridsize

nlat = toint((latN-latS)/dlat) + 1
mlon = toint((lonE-lonW)/dlon)

lat  = fspan(todouble(latS), todouble(latN), nlat)
lon  = fspan(todouble(lonW), todouble(lonE)-todouble(dlon), mlon)

lat@units = "degrees_north"
lon@units = "degrees_east"

count     = new( (/nlat,mlon/), "integer", -99)
count!0   = "lat"
count!1   = "lon"
count&lat =  lat
count&lon =  lon
count=0



;  shp_filename="../shapefile_ncl/statesp010g.shp"
;  opt = True
;  opt@minlat=37.
;  opt@maxlat=49.
;  opt@minlon=-82.
;  opt@maxlon=-66.
;  opt@shape_var = "NAME"
;  opt@DEBUG="True"

if (masktype .eq. "NCA") then
  sname = "shapefiles/NCA_polygon.shp"
  shapevarname = "NAME"
  seperate_regions=(/"Northeast","East North Central","West North Central","Northwest","Southwest","South","Southeast"/)
  ;seperate_regions=(/"East North Central","West North Central"/)
else if (masktype .eq. "RSI") then
  sname = "/Users/cmz5202/Software/esta/calc_RSI/shapefile_ncl/statesp010g.shp"
  shapevarname = "NAME"
  seperate_regions=(/"Northeast","Central","East North Central","West North Central","South","Southeast"/)
else
  print("none")
end if
end if

f = addfile(sname,"r")
print_shapefile_info(sname)
plot_shapefile(sname)

do ii = 0,dimsizes(seperate_regions)-1

  if (masktype .eq. "RSI") then
    if (seperate_regions(ii) .eq. "Northeast") then
      region_arr = (/"Pennsylvania","New Jersey","Delaware","Maryland","Massachusetts","New York","Connecticut","Vermont","Maine","Rhode Island","New Hampshire","District of Columbia"/)
    elseif (seperate_regions(ii) .eq. "Central") then
      region_arr = (/"Ohio","Indiana","Illinois","Missouri","Kentucky","Tennessee","West Virginia"/)
    elseif (seperate_regions(ii) .eq. "East North Central") then
      region_arr = (/"Michigan","Wisconsin","Iowa","Minnesota"/)
    elseif (seperate_regions(ii) .eq. "West North Central") then
      region_arr = (/"Montana","North Dakota","South Dakota","Wyoming","Nebraska"/)
    elseif (seperate_regions(ii) .eq. "South") then
      region_arr = (/"Kansas","Texas","Oklahoma","Arkansas","Louisiana","Mississippi"/)
    elseif (seperate_regions(ii) .eq. "Southeast") then
      region_arr = (/"Alabama","Georgia","Florida","South Carolina","North Carolina","Virginia"/)
    end if
  else
    region_arr = (/seperate_regions(ii)/)
  end if

  opt=True
  opt@shape_var=(/shapevarname/)
  opt@shape_names=region_arr
  opt@return_mask=True
  opt@debug=True

  maskdata = shapefile_mask_data(count,sname,opt)
  maskdata = maskdata*(ii+1)
  count = count+maskdata

  delete(opt)
  delete(region_arr)

end do


system("/bin/rm -f simple.nc")   ; remove any pre-existing file
ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Simple Approach to netCDF Creation"
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

ncdf->maskdata = count


;0       1           Northeast  MULTIPOLYGON (((-77.18389 38.57444, -77.18166 ...
;1       2  East North Central  POLYGON ((-94.61816 36.49837, -94.61847 36.600...
;2       3  West North Central  POLYGON ((-95.78650 40.00150, -95.99983 40.001...
;3       4           Northwest  MULTIPOLYGON (((-124.21111 43.40167, -124.2102...
;4       5           Southwest  MULTIPOLYGON (((-119.82500 33.95196, -119.8247...
;5       6               South  MULTIPOLYGON (((-95.07222 29.78861, -95.07166 ...
;6       7           Southeast  MULTIPOLYGON (((-81.92486 24.53308, -81.92486 ...
