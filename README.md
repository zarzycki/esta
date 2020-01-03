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

**A sample of gridded reanalysis data to start with can be downloaded from:** http://www.personal.psu.edu/cmz5202/data/esta-JRA-sample.tar.gz

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
SNOW (snowfall, averaged or accumulated over prior data step, arbitrary units)

Data should be in NetCDF format,on a regular latitude-longitude grid, and must have an associated time dimension.

PSL (time, lat, lon)
SNOW (time, lat, lon)

- There must be a coordinate variable lat monotonically increasing from S->N
- There must be a coordinate variable lon monotonically increasing from W->E
- Time must be in a compliant format (ex: "hours since 2000-01-01 00:00") and contain both a "units" and "calendar" attribute).
- Data should be 6-hourly (00Z, 06Z, 12Z, 18Z), although other increments are supported with minor tweaks to the code. However, PSL and SNOW data should be at the same temporal resolution (i.e., if 3-hourly PSL data is available but 6-hourly SNOW, PSL should be subsampled every other time to match).

Sample data files for January 2016 (JRA-55) can be downloaded:

```
wget http://www.personal.psu.edu/cmz5202/data/esta-JRA-sample.tar.gz
```

If you have a dataset you believe to be a standard, compliant format but breaks ESTA, please contact me.

### 2.) Build filelist lookups

TempestExtremes and ESTA use file lookup tables to track storms across multiple files when necessary. These tables essentially concatenate along the record (time) dimension.

Two files should be generated, with full path to all files to be analyzed store one per line. One file should contain the series of files where the code can find PSL. One file should contain the series of files where the code can find SNOW. Note that files *may* contain both PSL and SNOW; in this case only one filelist is necessary and should be appropriately referenced in the namelist.

An example script to generate a filelist is: $ESTAPATH

### 3.) Set namelist variables

The namelist is the "conductor" for ESTA. Many options need to be set, but the hope is by setting them correctly, the code itself can function without being modified by the end-user. The table of options is below. A sample namelist which can be copied and then modified for new data can be found at XXXXX.

| Namelist Variable | Namelist sample | Type | Description |
| --- | --- | --- | --- |
| ESTAPATH | ESTAPATH="/home/zarzycki/esta/" | string | Path to main ESTA dir |
| TE_TEMPESTEXTREMESDIR | TE_TEMPESTEXTREMESDIR="/home/zarzycki/TE/" | string | Path to main TempestExtremes dir |
| TE_UQSTR | TE_UQSTR="JRA.2016" | string | Unique string for TE traj calcs |
| TE_CONNECTFLAG | TE_CONNECTFLAG="" | string | TE connectivity string for unstructured grids (empty if unused) |
| TE_MPIRUNCMD | TE_MPIRUNCMD="" | string | TE MPI run command (empty if no MPI) |
| FILELISTNAME | FILELISTNAME="" | string | List of files containing PSL fields for ETC tracking |
| TRAJDIR | TRAJDIR="./" | string | Directory to output ETC trajectory file |
| LIST_SNOW | LIST_SNOW="filelist.SNOW.JRA.2016","filelist.PSL.JRA.2016" | string(s) | List of files containing event variables to be extracted |
| DELTADEG | DELTADEG=12.0 | float | Radius of precipitation integral search (in great circle degrees) |
| SWE | SWE="12" | string | Snow water equivalent for snowfall calculations |
| EX_OUTDIR | EX_OUTDIR="./storms/" | string | Path to output storm extraction NetCDF file |
| EX_OUTFILE | EX_OUTFILE="JRA.2016.nc" | string | Name of storm extraction NetCDF file |
| EX_INVARS | EX_INVARS="SNOW","PSL" | string(s) | Name of variables to extract |
| EX_AGGOP | EX_AGGOP="sum","min" | string(s) | Integration operator for storm extraction |
| EX_OFFSET | EX_OFFSET=1,1 | float(s) | Multiplication operators for extract variables |
| EX_DOTIMESERIES | EX_DOTIMESERIES=True,True | logical(s) | Multiplication operators for extract variables |
| EX_DOREGOUT | EX_DOREGOUT=True,True | logical(s) | Only output regional NetCDF domain for storm extraction |
| EX_FORCEDELOUT | EX_FORCEDELOUT=True,False | logical(s) | If EX_OUTFILE exists, should we force delete? |
| RSI_SUMVAR | RSI_SUMVAR="SUM_SNOW" | string | Name of variable to sum over for RSI |
| RSI_OUTDIR | RSI_OUTDIR="./RSI/" | string | Where to write RSI output |

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



