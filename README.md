# ESTA pre-tag

## a.k.a. the Extratropical Snowstorm Tracking Algorithm

Colin M. Zarzycki

ESTA allows for extratropical cyclones (ETCs) to be tracked in gridded weather and climate data. Snowfall (or other precipitation) files are extracted following the ETC center and integrated to produce storm totals. An NCL tool allows for storm-integrated snowfall to be overlaid with a population density map to produce Regional Snowfall Index (RSI) projections.

This code is used in the paper:
_Zarzycki, C. M. (2018). Projecting changes in societally impactful northeastern U.S. snowstorms. Geophysical Research Letters, 45, 12,067–12,075. https://doi.org/10.1029/2018GL079820._

**NOTE:** This is pre-release code that mainly consists of NCL / Bash scripting.

**WARNING:** This code has not been extensively verified beyond the particular project noted above. If you find any bugs or inconsistencies, please contact czarzycki [at] psu [dot] edu.

It applies the **Regional Snowfall Index** as defined by:
_Squires, M.F., J.H. Lawrimore, R.R. Heim, D.A. Robinson, M.R. Gerbush, and T.W. Estilow, 2014: The Regional Snowfall Index. Bull. Amer. Meteor. Soc., 95, 1835–1848, https://doi.org/10.1175/BAMS-D-13-00101.1._

To use this code, you *must* have gridded reanalysis or climate model data that contains at least sea level pressure and snowfall (averaged over the previous "timestep" of your data. This is preferred and well-supported at 6-hourly and finer, but 12-hourly and 24-hourly may still be of use, particularly for more intense ETCs.

**A sample of gridded reanalysis data to start with can be downloaded from:** ?????

## General procedure

1. Pre-process sea level pressure and precipitation data to appropriate format and generate file/path list in txt format.
2. Set namelist variables.
3. Build filelist lookups.
4. Run tracker (find ETCs).
5. Run extractor (pull out storm-integrated quantities).
6. Run RSI analysis (overlay storm-integrated snowfall w/ population density).

## Detailed procedure

Notes:
- ESTAPATH is the 
- sdf


### 1.) Pre-process data

Two variables are needed for ESTA to properly track storms:
PSL (sea level pressure, in Pa)
SNOW (snowfall, averaged over prior data step, arbitrary units)

Data should be in NetCDF format,on a regular latitude-longitude grid, and must have an associated time dimension.

PSL (time, lat, lon)
SNOW (time, lat, lon)

- There must be a coordinate variable lat monotonically increasing from S->N
- There must be a coordinate variable lon monotonically increasing from W->E
- Time must be in a compliant format (ex: "hours since 2000-01-01 00:00") and contain both a "units" and "calendar" attribute).
- Data should be 6-hourly (00Z, 06Z, 12Z, 18Z), although other increments are supported with minor tweaks to the code. However, PSL and SNOW data should be at the same temporal resolution (i.e., if 3-hourly PSL data is available but 6-hourly SNOW, PSL should be subsampled every other time to match).

See example of data files in XXXXXX. If you have a dataset you believe to be a standard, compliant format but breaks ESTA, please contact me.

### 2.) Build filelist lookups

TempestExtremes and ESTA use file lookup tables to track storms across multiple files when necessary. These tables essentially concatenate along the record (time) dimension.

Two files should be generated, with full path to all files to be analyzed store one per line. One file should contain the series of files where the code can find PSL. One file should contain the series of files where the code can find SNOW. Note that files *may* contain both PSL and SNOW; in this case only one filelist is necessary and should be appropriately referenced in the namelist.

An example script to generate a filelist is: $ESTAPATH

### 3.) Set namelist variables

The namelist is the "conductor" for ESTA. Many options need to be set, but the hope is by setting them correctly, the code itself can function without being modified by the end-user. The table of options is below. A sample namelist which can be copied and then modified for new data can be found at XXXXX.

### 4.) Run tracker

```
${ESTAPATH}/tracking/esta-etc-tracking.sh ${NAMELIST} 
```

### 5.) Run extractor

```
${ESTAPATH}/tracking/esta-extract-storms.sh ${NAMELIST} 
```

### 6.) Run RSI analysis

```
${ESTAPATH}/tracking/esta-calc-RSI.sh ${NAMELIST}
```



