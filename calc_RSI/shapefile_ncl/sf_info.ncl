;----------------------------------------------------------------------
; This script prints basic information about a shapefile. To execute:
;
;   ncl 'sname="xxxx.shp"' sf_info.ncl
;
; where "xxxx.shp" is the name of the shapefile
;
;----------------------------------------------------------------------

begin
  dq = str_get_dq()                   ; double quote 

  if(.not.isvar("sname")) then
    print("To run this script, you must provide a filename:")
    print("   ncl 'sname=" + dq + "xxxx.shp" + dq + "' sf_info.ncl")
    exit
  end if

;---Open shapefile
  f = addfile(sname,"r")

;---Read data off shapefile
  var_names    = getfilevarnames(f)                  ; variable names
  geo_dims     = getfilevardimsizes(f,"geometry")    ; features on file
  num_features = geo_dims(0)
  lon          = f->x
  lat          = f->y
  nvars        = dimsizes(var_names)

  print("======================================================================")
  print("Filename: " + dq + sname + dq)
  print("   Geometry type: " + f@geometry_type)
  print("   # of features: " + num_features) 
  print("   Min/max lat:   " + sprintf("%7.2f",min(lat)) + "/" + sprintf("%7.2f",max(lat)))
  print("   Min/max lon:   " + sprintf("%7.2f",min(lon)) + "/" + sprintf("%7.2f",max(lon)))
  print("   Variable names and their types:")
  do n=0,nvars-1
     print("       " + var_names(n) + " : " + getfilevartypes(f,var_names(n)))
  end do
  print("======================================================================")
end
