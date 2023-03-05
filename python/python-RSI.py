import geopandas 
import xarray
import rasterio.features
import fiona
import numpy as np

shapefile = './shapefiles/NCA_polygon.shp'

sf = geopandas.read_file(shapefile)
print(sf)



US=sf[sf['NAME']=='Northeast']  # conus, AK, HI, PR

print(US)


coords = US['geometry']['x'][0]
print(coords)

