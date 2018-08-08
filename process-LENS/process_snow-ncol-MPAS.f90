! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This is an example which reads some 4D pressure and
! temperatures. The data file read by this program is produced by
! the companion program pres_temp_4D_wr.f90. It is intended to
! illustrate the use of the netCDF Fortran 90 API.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: pres_temp_4D_rd.f90,v 1.6 2006/12/09 18:44:58 russ Exp $

program pres_temp_4D_rd
  use netcdf
  implicit none

! add an interface block for DRELHUM ...
	interface
		REAL FUNCTION DRELHUM(T,W,P)
      REAL T
      REAL W
      REAL P
		end function DRELHUM
	end interface


  CHARACTER(len=256) :: unitstring, timeperstring, indir, outdir
  CHARACTER(len=1) :: ensmem

  ! This is the name of the data file we will read.
  character (len = 256) :: TFILENAME, UFILENAME, VFILENAME, QFILENAME, ZFILENAME
  character (len = 256) :: PSFILENAME, PRECTFILENAME

  integer :: t_id, u_id, v_id, q_id, z_id, prect_id, ps_id, pm_id, pi_id
  integer :: ncid

  ! We are reading 4D data, a 2 x 6 x 12 lev-lat-lon grid, with 2
  ! timesteps of data.

!! NE30 = 48602
!! NE0CONUS = 174098
!! MP120A = 40962
!! MP15A-120A = 195266
  integer, parameter :: NDIMS = 2, NRECS = 120
  integer, parameter :: NLEVS = 32, NLATS = 195266
  character (len = *), parameter :: LEV_NAME = "lev"
  character (len = *), parameter :: LAT_NAME = "ncol"
  character (len = *), parameter :: REC_NAME = "time"
  integer :: lev_dimid, lat_dimid, rec_dimid

  ! Constants
  real, parameter :: P0 = 100000.

  ! The start and count arrays will tell the netCDF library where to
  ! read our data.
  integer :: start(NDIMS), count(NDIMS)

  ! In addition to the latitude and longitude dimensions, we will also
  ! create latitude and longitude variables which will hold the actual
  ! latitudes and longitudes. Since they hold data about the
  ! coordinate system, the netCDF term for these is: "coordinate
  ! variables."
  real :: lats(NLATS), levs(NLEVS), recs(NRECS)
  integer :: lat_varid, lev_varid, rec_varid

  ! Loop integers
  integer :: i,j,k,z
  integer :: IER, IOPT

  ! variables for playing with weightings of profiles +/- 6 hours.
  integer :: HOPT
  real :: W1, W2

  ! random seed for random number gen
  integer :: rdseed

  ! We will read surface temperature and pressure fields. In netCDF
  ! terminology these are called "variables."
  integer :: hyai_varid,hybi_varid,hyam_varid,hybm_varid
  integer :: prect_varid,ps_varid,temp_varid,u_varid,v_varid,q_varid,z3_varid,pm_varid,pi_varid
  integer :: prect_snow_varid,prect_ice_varid,prect_fzra_varid
  integer :: prect_mix_varid,prect_rain_varid
  integer :: vort_varid, div_varid, omega_varid, rh_varid
  integer :: dimids(NDIMS),dimids4d(3)
  integer :: ptype_varid,ptypecz_varid,ratio_varid,snowrate_varid

  integer :: pmid_varid, pdel_varid, psx_varid, psy_varid, psxin_varid, psyin_varid

  ! Program variables to hold the data we will read in. We will only
  ! need enough space to hold one timestep of data; one record.
  real :: hyai(NLEVS+1)
  real :: hybi(NLEVS+1)
  real :: hyam(NLEVS)
  real :: hybm(NLEVS)
  real :: temp_in(NLATS,NLEVS,NRECS)
  real :: prect_in(NLATS,NRECS)
  real :: ps_in(NLATS,NRECS)
  !real :: u(NLATS,NLEVS,NRECS)
  !real :: v(NLATS,NLEVS,NRECS)
  !real :: q(NLATS,NLEVS,NRECS)
  real :: z3(NLATS,NLEVS,NRECS)
  
  real :: prect_snow(NLATS,NRECS)
  real :: prect_rain(NLATS,NRECS)
  real :: prect_ice(NLATS,NRECS)
  real :: prect_fzra(NLATS,NRECS)
  real :: prect_mix(NLATS,NRECS)
  real :: psgx(NLATS,NRECS)
  real :: psgy(NLATS,NRECS)
  !real :: vort(NLATS,NLEVS,NRECS)
  !real :: div(NLATS,NLEVS,NRECS)
  !real :: omega(NLATS,NLEVS,NRECS)
  !real :: relhum(NLATS,NLEVS,NRECS)

  !real :: snowratio(NLATS,NRECS)
  real :: snowfallrate(NLATS,NRECS)

  real :: pdel(NLATS,NLEVS,NRECS)
  real :: pmid(NLATS,NLEVS,NRECS)
  real :: pint(NLATS,NLEVS+1,NRECS)
  real :: zint(NLATS,NLEVS+1,NRECS)

  ! Variables we modify here
  real :: pres_m(NLEVS)
  real :: pres_i(NLEVS+1)
  real :: hybd(NLEVS)
  real :: z_i(NLEVS+1)

  real :: logps(NLATS,NRECS)

  integer :: ptype(NLATS,NRECS)
  integer :: ptypecz(NLATS,NRECS)

  integer :: nprlev
  real, dimension(NLEVS) :: dumpmid, dumz3

  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = 128) :: LAT_UNITS, LEV_UNITS, REC_UNITS, PRECT_UNITS

  ! To check the units attributes.
  character*80 pres_units_in, temp_units_in
  character*80 lat_units_in

  character (len = 256) :: FILE_NAME

  integer :: stat

  ! Get command line arguments and set filenames
  CALL getarg(1, indir)
  CALL getarg(2, outdir)
  CALL getarg(3, ensmem)
  !print *,'ensmem: ',ensmem,'unitstring: ',unitstring
  print *,'indir: ',indir
  print *,'outdir: ',outdir
  print *,'option: ',ensmem
  read(ensmem,*,IOSTAT=stat) HOPT
  print *, HOPT

  TFILENAME = trim(""//trim(indir)//"")
  print *,TFILENAME

!   UFILENAME = trim(""//trim(indir)//"")
! 
!   VFILENAME = trim(""//trim(indir)//"")
! 
!   QFILENAME = trim(""//trim(indir)//"")

  ZFILENAME = trim(""//trim(indir)//"")

  PSFILENAME = trim(""//trim(indir)//"")

  PRECTFILENAME = trim(""//trim(indir)//"")

  FILE_NAME = trim(""//trim(outdir)//"")

  ! Open the file. 
  call check( nf90_open(TFILENAME, nf90_nowrite, t_id) )
!   call check( nf90_open(QFILENAME, nf90_nowrite, q_id) )
!   call check( nf90_open(UFILENAME, nf90_nowrite, u_id) )
!   call check( nf90_open(VFILENAME, nf90_nowrite, v_id) )
  call check( nf90_open(ZFILENAME, nf90_nowrite, z_id) )
  call check( nf90_open(PRECTFILENAME, nf90_nowrite, prect_id) )
  call check( nf90_open(PSFILENAME, nf90_nowrite, ps_id) )

  print *,'1'

  ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(t_id, "lev", lev_varid) )
  call check( nf90_inq_varid(t_id, "time",rec_varid) )

  print *,'2'

  ! Read the latitude and longitude data.
  call check( nf90_get_var(t_id, lev_varid, levs) )
  call check( nf90_get_var(t_id, rec_varid, recs) )

  print *,'3'

  call check( nf90_get_att(t_id, lev_varid, UNITS, LEV_UNITS) )
  call check( nf90_get_att(t_id, rec_varid, UNITS, REC_UNITS) )

  print *,'4'

  print *,'5'

  call check( nf90_inq_varid(prect_id, "PRECT", prect_varid  ))
  call check( nf90_inq_varid(ps_id, "PS",    ps_varid ))
  call check( nf90_inq_varid(t_id, "T",     temp_varid) )
!   call check( nf90_inq_varid(u_id, "U",     u_varid) )
!   call check( nf90_inq_varid(v_id, "V",     v_varid) )
!   call check( nf90_inq_varid(q_id, "Q",     q_varid) )
  call check( nf90_inq_varid(z_id, "Z3",     z3_varid) )
  call check( nf90_inq_varid(t_id, "PRESSURE",     pm_varid) )
  call check( nf90_inq_varid(t_id, "PRESSUREi",    pi_varid) )

  print *,'6'

  call check( nf90_get_att(prect_id, prect_varid, UNITS, PRECT_UNITS) )

  print *,'7'

  call check( nf90_get_var(prect_id, prect_varid, prect_in ) )
  call check( nf90_get_var(ps_id, ps_varid,    ps_in ) )
  call check( nf90_get_var(t_id, temp_varid,  temp_in) )
  !call check( nf90_get_var(u_id, u_varid,     u) )
  !call check( nf90_get_var(v_id, v_varid,     v) )
  !call check( nf90_get_var(q_id, q_varid,     q) )
  call check( nf90_get_var(z_id, z3_varid,     z3) )
  call check( nf90_get_var(t_id, pm_varid,     pmid) )
  call check( nf90_get_var(t_id, pi_varid,     pint) )

  print *,'8'

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(t_id) )
!   call check( nf90_close(u_id) )
!   call check( nf90_close(v_id) )
!   call check( nf90_close(q_id) )
  call check( nf90_close(z_id) )
  call check( nf90_close(prect_id) )
  call check( nf90_close(ps_id) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS reading files !"

  !prect_in = prect_in*141732.283
  !call test(temp_in)

!!! DO STUFF HERE

  ptype = 0
  ptypecz = 0
  !omega = 0
  pdel = 0
  logps = LOG(ps_in)
  !call bourgouin(hyam,hyai,hyam,hybi,ptype)


! (NLONS,NLATS,NLEVS,NRECS)
!  do z = 1,NRECS
!    !do k = 1, NLVLS
!      do i = 1, NLATS
!        do j = 1, NLONS         
!          pres_m = hyam*P0+hybm*ps_in(j,i,z)
!          pres_i = hyai*P0+hybi*ps_in(j,i,z)
!          call bourgouin(pres_m,pres_i,temp_in(j,i,:,z),ptype)
!        end do
!      end do
!    !end do
!  end do

print *,'Getting interfaces and midpoints...'
  do z = 1,NRECS
      do i = 1,NLATS
          pres_m = pmid(i,:,z)
          pres_i = pint(i,:,z)
          pdel(i,:,z) = pres_i(2:size(pres_i))-pres_i(1:size(pres_i)-1)

          dumpmid = pmid(i,:,z)
          dumz3 = z3(i,:,z)
          call DINT2P(pres_m,z3(i,:,z),dumpmid,dumz3,NLEVS,pres_i,z_i,NLEVS+1,-1,1.0E+36,IER)
          zint(i,:,z) = z_i
      end do
  end do
print *,'... interfaces/midpoints calc done'


print *,'Starting ptype calculations...'
  do j = 1,NRECS
     !call grad2DCFD(logps(:,:,z),psgx(:,:,z),psgy(:,:,z),lats,lons)
     !call DPHYBRID(P0,hyam,hybm,ps_in(:,:,z),NLONS,NLATS,NLEVS,pmid(:,:,:,z),IER)
     !call DDPHYBRID(P0,hyai,hybi,ps_in(:,:,z),NLONS,NLATS,NLEVS,pdel(:,:,:,z),IER)
      do i = 1,NLATS
        !IF (mod(i,100) .EQ. 0) print *,i,' of ',NLATS
          z=NRECS-j+1
          rdseed = i*z
            IF (HOPT .GT. 0) THEN
              IF (HOPT .EQ. 1) THEN     
                W1 = 1.
                W2 = 1.
              ELSE IF (HOPT .EQ. 2) THEN
                W1 = 2.
                W2 = 1.
              ELSE IF (HOPT .EQ. 3) THEN
                W1 = 1.
                W2 = 2.
              END IF 
              IF (j .EQ. 1 .AND. i .EQ. 1) THEN
                print *,'Doing weighting!'
                print *,'W1: ',W1
                print *,'W2: ',W2
              END IF
              IF (z .NE. 1) THEN
                temp_in(i,:,z)=(W1*temp_in(i,:,z)+W2*temp_in(i,:,z-1))/(W1+W2)
                pmid(i,:,z)=(W1*pmid(i,:,z)+W2*pmid(i,:,z-1))/(W1+W2)
                pint(i,:,z)=(W1*pint(i,:,z)+W2*pint(i,:,z-1))/(W1+W2)
                zint(i,:,z)=(W1*zint(i,:,z)+W2*zint(i,:,z-1))/(W1+W2)
              END IF
            END IF

          call bourgouin(rdseed,pmid(i,:,z),pint(i,:,z),temp_in(i,:,z),ptypecz(i,z))
          call calwxt_bourg(rdseed,32,33,temp_in(i,:,z),pmid(i,:,z),pint(i,:,z),zint(i,:,z),ptype(i,z))
          
          !do k = 1,NLEVS
          !  relhum(i,k,z)=DRELHUM(temp_in(i,k,z),q(i,k,z),pmid(i,k,z))          
          !end do
      end do
  end do
  
  nprlev = 0
  !do k=1,NLEVS
  !  if (nprlev .EQ. 0 .and. hybi(k) .GT. 0.0) then
  !    nprlev = k
  !  end if
  !end do
  !nprlev = 12
print *,'... ptype calcs complete'
  

  ! init arrays to zero
  !vort = 0
  !div = 0
  !snowratio = 0

!>   ! Calculate vorticity and divergence
!> print *,'Starting VORT/DIV calc...'
!>   do z = 1,NRECS
!>   do k = 1,NLEVS
!>     call DVRFIDF(u(:,:,k,z),v(:,:,k,z),lats,lons,NLONS,NLATS,0,IOPT,vort(:,:,k,z),IER)
!>     call DDVFIDF(u(:,:,k,z),v(:,:,k,z),lats,lons,NLONS,NLATS,0,IOPT, div(:,:,k,z),IER)
!>     IF (IER .NE. 0) print *,'IER NE ZERO at k=',k,' z=',z
!>   end do
!>   end do
!> print *,'... VORT/DIV calc done'
!> 
!> print *,'Starting OMEGA calc...'
!>   do z = 1,NRECS
!>     call omcalcccm(u(:,:,:,z),v(:,:,:,z),div(:,:,:,z),psgx(:,:,z),psgy(:,:,z),pmid(:,:,:,z),pdel(:,:,:,z),&
!>       ps_in(:,:,z),hybd,hybm,nprlev,omega(:,:,:,z),NLONS,NLATS,NLEVS)
!>   end do
!> print *,'... OMEGA calc done'
!> 
!> print *,'Starting Cobb calc...'
!>   do z = 1,NRECS
!>       do i = 1,NLATS
!>         do j = 1,NLONS
!>           call cobb(omega(j,i,:,z),temp_in(j,i,:,z),relhum(j,i,:,z),zint(j,i,:,z),pmid(j,i,:,z),snowratio(j,i,z))
!>         end do
!>       end do
!>   end do
!> print *,'... Cobb calc done'

print *,'Starting mask and rate calculations...'

  ! mask off the different ptypes
  prect_snow = prect_in
  call maskprecip3D(ptype,prect_snow,0)
  prect_mix = prect_in
  call maskprecip3D(ptype,prect_mix,1)
  prect_rain = prect_in
  call maskprecip3D(ptype,prect_rain,2)
  prect_ice = prect_in
  call maskprecip3D(ptype,prect_ice,3)
  prect_fzra = prect_in
  call maskprecip3D(ptype,prect_fzra,4)

print *,'... mask/rate calcs done'


print *,'OUTPUT FILE'

  call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )

  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  call check( nf90_def_dim(ncid, "lev", NLEVS, lev_dimid) )
  call check( nf90_def_dim(ncid, "ncol", NLATS, lat_dimid) )
  call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, rec_dimid) )

  ! Define the coordinate variables. We will only define coordinate
  ! variables for lat and lon.  Ordinarily we would need to provide
  ! an array of dimension IDs for each variable's dimensions, but
  ! since coordinate variables only have one dimension, we can
  ! simply provide the address of that dimension ID (lat_dimid) and
  ! similarly for (lon_dimid).
  call check( nf90_def_var(ncid, "ncol", NF90_REAL, lat_dimid, lat_varid) )
  call check( nf90_def_var(ncid, "lev", NF90_REAL, lev_dimid, lev_varid) )
  call check( nf90_def_var(ncid, "time", NF90_REAL, rec_dimid, rec_varid) )

  ! Assign units attributes to coordinate variables.
  call check( nf90_put_att(ncid, lat_varid, UNITS, trim(LAT_UNITS)) )
  call check( nf90_put_att(ncid, lev_varid, UNITS, trim(LEV_UNITS)) )
  call check( nf90_put_att(ncid, rec_varid, UNITS, trim(REC_UNITS)) )
  call check( nf90_put_att(ncid, rec_varid, "calendar", "noleap") )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids = (/ lat_dimid, rec_dimid /)
  dimids4d = (/ lat_dimid, lev_dimid, rec_dimid /)

  ! Define the netCDF variables for the pressure and temperature data.
  call check( nf90_def_var(ncid, "PRECT", NF90_REAL, dimids, prect_varid) )
  call check( nf90_def_var(ncid, "PRECBSN", NF90_REAL, dimids, prect_snow_varid) )
  call check( nf90_def_var(ncid, "PRECBMX", NF90_REAL, dimids, prect_mix_varid) )
  call check( nf90_def_var(ncid, "PRECBRA", NF90_REAL, dimids, prect_rain_varid) )
  call check( nf90_def_var(ncid, "PRECBIP", NF90_REAL, dimids, prect_ice_varid) )
  call check( nf90_def_var(ncid, "PRECBFZ", NF90_REAL, dimids, prect_fzra_varid) )
!  call check( nf90_def_var(ncid, "VORT", NF90_REAL, dimids4d, vort_varid) )
!  call check( nf90_def_var(ncid, "DIV", NF90_REAL, dimids4d, div_varid) )
!  call check( nf90_def_var(ncid, "OMEGA", NF90_REAL, dimids4d, omega_varid) )
!  call check( nf90_def_var(ncid, "RH", NF90_REAL, dimids4d, rh_varid) )

! Debugs
!  call check( nf90_def_var(ncid, "pmid", NF90_REAL, dimids4d, pmid_varid) )
!  call check( nf90_def_var(ncid, "PS", NF90_REAL, dimids, ps_varid) )
!  call check( nf90_def_var(ncid, "pdel", NF90_REAL, dimids4d, pdel_varid) )
!  call check( nf90_def_var(ncid, "PSX", NF90_REAL, dimids, psx_varid) )
!  call check( nf90_def_var(ncid, "PSY", NF90_REAL, dimids, psy_varid) )

   
!  call check( nf90_def_var(ncid, "U850",  NF90_REAL, dimids, u850_varid) )
!  call check( nf90_def_var(ncid, "T",     NF90_REAL, dimids4d, temp_varid) )
  call check( nf90_def_var(ncid, "PTYPE",     NF90_REAL, dimids, ptype_varid) )
  call check( nf90_def_var(ncid, "PTYPECZ",     NF90_REAL, dimids, ptypecz_varid) )

  ! Assign units attributes to the netCDF variables.
  call check( nf90_put_att(ncid, prect_varid, UNITS, trim(PRECT_UNITS)) )
  call check( nf90_put_att(ncid, prect_snow_varid, UNITS, trim(PRECT_UNITS)) )
  call check( nf90_put_att(ncid, prect_mix_varid, UNITS, trim(PRECT_UNITS)) )
  call check( nf90_put_att(ncid, prect_rain_varid, UNITS, trim(PRECT_UNITS)) )
  call check( nf90_put_att(ncid, prect_ice_varid, UNITS, trim(PRECT_UNITS)) )
  call check( nf90_put_att(ncid, prect_fzra_varid, UNITS, trim(PRECT_UNITS)) )
  
  ! End define mode.
  call check( nf90_enddef(ncid) )
  
  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  call check( nf90_put_var(ncid, lat_varid, lats) )
  call check( nf90_put_var(ncid, lev_varid, levs) )
  call check( nf90_put_var(ncid, rec_varid, recs) )

  call check( nf90_put_var(ncid, prect_varid, prect_in))
  call check( nf90_put_var(ncid, prect_snow_varid, prect_snow))
  call check( nf90_put_var(ncid, prect_mix_varid, prect_mix))
  call check( nf90_put_var(ncid, prect_rain_varid, prect_rain))
  call check( nf90_put_var(ncid, prect_ice_varid, prect_ice))
  call check( nf90_put_var(ncid, prect_fzra_varid, prect_fzra))
  call check( nf90_put_var(ncid, ptype_varid, ptype))
  call check( nf90_put_var(ncid, ptypecz_varid, ptypecz))



!  call check( nf90_put_var(ncid, vort_varid, vort))
!  call check( nf90_put_var(ncid, div_varid, div))
!  call check( nf90_put_var(ncid, omega_varid, omega))
!  call check( nf90_put_var(ncid, rh_varid, relhum))

!  call check( nf90_put_var(ncid, pmid_varid, pmid))
!  call check( nf90_put_var(ncid, pdel_varid, pdel))
!  call check( nf90_put_var(ncid, psx_varid, psgx))
!  call check( nf90_put_var(ncid, psy_varid, psgy))
!  call check( nf90_put_var(ncid, ps_varid, ps_in))
  
  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are really written to disk.
  call check( nf90_close(ncid) )
  
  print *,"*** SUCCESS writing example file !"

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

subroutine test(var)
  real, intent(inout) :: var(NLATS,NLEVS,NRECS)
  var = var-273.15
end subroutine test

subroutine bourgouin(rdseed,pres_m,pres_i,T_m,ptype)

  implicit none

  real, intent(in) :: pres_m(NLEVS)
  real, intent(in) :: pres_i(NLEVS+1)
  real, intent(in) :: T_m(NLEVS)
  integer, intent(in) :: rdseed
  integer, intent(inout) :: ptype

  ! Declare local variables
  real :: pres_i_top(NLEVS)
  real :: pres_i_bot(NLEVS)
  real :: area(NLEVS)
  integer :: area_sign(NLEVS)
  integer :: signchange
  integer :: counter
  integer :: rseed
  real :: r1  ! random number
  REAL, DIMENSION(:), ALLOCATABLE :: area_arr,dumareaPA,dumareaNA
  real :: PA_a, PA_sfc, NA, PA

  ! Loop indices
  integer :: i,j,k
  integer :: AllocateStatus,DeAllocateStatus

  ! Declare local constants
  REAL, PARAMETER :: eps = 10E-10
  REAL, PARAMETER :: Rd = 287.0
  REAL, PARAMETER :: cp = 1004.0
 
  ! set random seed
  rseed=rdseed
  CALL RANDOM_SEED(rseed)

  pres_i_top = pres_i(1:NLEVS)
  pres_i_bot = pres_i(2:NLEVS+1)

  !print *,pres_i_top


  area = -cp*(T_m-273.15)*log(pres_i_top/pres_i_bot)
!  area = where(abs(area).le.eps,eps,area)

  WHERE(abs(area)<eps)
    area = eps
  ELSEWHERE
    area = area 
  END WHERE

  !print *,area

  area_sign = area/abs(area)

  ! Find number of sign changes
  signchange = 0
  DO k = 1,SIZE(area_sign)-1
    IF (area_sign(k) - area_sign(k+1) .NE. 0) THEN
      signchange = signchange+1
    END IF
  END DO


  IF (signchange .eq. 0) THEN
    ! SNOW
    ptype = 0
  ELSE
    !== We need to calculate PA/NA
    ALLOCATE( area_arr(signchange+1), STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
    area_arr=0
    counter=1
    DO k = 1,SIZE(area_sign)-1
      area_arr(counter) = area_arr(counter)+area(k)
      IF (area_sign(k) - area_sign(k+1) .NE. 0) THEN
        counter = counter+1
      END IF
    END DO
    ! Need to add last bit?
    area_arr(counter) = area_arr(counter)+area(SIZE(area_sign))

    IF (signchange .eq. 1) THEN
      PA = area_arr(2)
      IF (PA.lt.5.6) THEN
        ! === SNOW
        ptype = 0
      ELSE IF (PA.gt.13.2) THEN
        ! === RAIN
        ptype = 2
      ELSE
        !! === MIX
        !ptype = 1
        CALL RANDOM_NUMBER(r1)
        IF (r1.gt.0.5) THEN
          ptype = 0
        ELSE
          ptype = 2
        END IF
      END IF
    ELSE IF (signchange .eq. 2) THEN
      ! area_arr(0) wasted stuff above first cross
      ! area_arr(1) PA
      ! area_arr(2) NA
      PA = area_arr(2)
      NA = -area_arr(3)
      !print("PA: "+PA+"    NA: "+NA)
      IF (NA .gt. 66+0.66*PA) THEN
        ! ==== ICE
        ptype = 3
      ELSE IF (NA .lt. 46+0.66*PA) THEN
        ! ==== FZRA
        ptype = 4
      ELSE
        !! ==== MIX OF FZRA/ICE
        !ptype = 1
        CALL RANDOM_NUMBER(r1)
        IF (r1.gt.0.5) THEN
          ptype = 3
        ELSE
          ptype = 4
        END IF
      END IF
    ELSE IF (signchange .eq. 3) THEN
      ! area_arr(0) wasted stuff above first cross
      ! area_arr(1) PA_a
      ! area_arr(2) NA
      ! area_arr(3) PA_sfc
      PA_a = area_arr(2)
      NA = -area_arr(3)
      PA_sfc = area_arr(4)
      if (PA_a .gt. 2) then
      ! we use Eqn 4
        if (NA .gt. 66+0.66*PA_sfc) then
          ! ==== ICE
          ptype = 3
        else if (NA .lt. 46+0.66*PA_sfc) then
          ! ==== RAIN
          ptype = 2
        else
          !! ==== MIX OF ICE/RAIN
          !ptype = 1
          CALL RANDOM_NUMBER(r1)
          IF (r1.gt.0.5) THEN
            ptype = 3
          ELSE
            ptype = 2
          END IF
        end if
      else
        if (PA_sfc.lt.5.6) then
          ! === ICE
          ptype = 3
        else if (PA_sfc.gt.13.2) then
          ! === RAIN
          ptype = 2
        else
          !! === MIX
          !ptype = 1
          CALL RANDOM_NUMBER(r1)
          IF (r1.gt.0.5) THEN
            ptype = 3
          ELSE
            ptype = 2
          END IF
        end if
      end if
    else
      if (mod(signchange,2) .eq. 0) then
        area_arr(1)=0
      
        ! Need to allocate dummy arrays
        ALLOCATE( dumareaPA(signchange+1), STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        ALLOCATE( dumareaNA(signchange+1), STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        dumareaPA = area_arr
        dumareaNA = area_arr
        WHERE(area_arr .LT. 0) dumareaPA = 0
        WHERE(area_arr .GT. 0) dumareaNA = 0
        PA = sum(dumareaPA)
        NA = -sum(dumareaNA)

          if (NA .gt. 66+0.66*PA) then
            !; ==== ICE
            ptype = 3
            !print("ICE")
          else if (NA .lt. 46+0.66*PA) then
            !; ==== FZRA
            ptype = 4
            !print("FRZRA")
          else
            !; ==== MIX OF FZRA/ICE
            !ptype = 1
          CALL RANDOM_NUMBER(r1)
          IF (r1.gt.0.5) THEN
            ptype = 3
          ELSE
            ptype = 4
          END IF
          end if
      else 
        area_arr(1)=0
        PA_sfc = area_arr(SIZE(area_arr))
        area_arr(SIZE(area_arr)) = 0

        ! Need to allocate dummy arrays
        ALLOCATE( dumareaPA(signchange+1), STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        ALLOCATE( dumareaNA(signchange+1), STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        dumareaPA = area_arr
        dumareaNA = area_arr
        WHERE(area_arr .LT. 0) dumareaPA = 0
        WHERE(area_arr .GT. 0) dumareaNA = 0
        PA_a = sum(dumareaPA)
        NA = -sum(dumareaNA)

        !;print("PA_a: "+PA_a+"    NA: "+NA+"    PA_sfc: "+PA_sfc)
        !;print(area)
        if (PA_a .gt. 2) then
        !; we use Eqn 4
          if (NA .gt. 66+0.66*PA_sfc) then
            !; ==== ICE
            ptype = 3
          else if (NA .lt. 46+0.66*PA_sfc) then
            !; ==== RAIN
            ptype = 2
          else
            !!; ==== MIX OF ICE/RAIN
            !ptype = 1
            CALL RANDOM_NUMBER(r1)
            IF (r1.gt.0.5) THEN
              ptype = 3
            ELSE
              ptype = 2
            END IF
          end if
        else
          if (PA_sfc.lt.5.6) then
            !; === ICE
            ptype = 3
          else if (PA_sfc.gt.13.2) then
            !; === RAIN
            ptype = 2
          else
            !; === MIX
            !ptype = 1
            CALL RANDOM_NUMBER(r1)
            IF (r1.gt.0.5) THEN
              ptype = 3
            ELSE
              ptype = 2
            END IF
          end if
        end if
        DEALLOCATE (dumareaPA, STAT = DeAllocateStatus)
        IF (DeAllocateStatus /= 0) STOP "*** Deallocate failed ***"
        DEALLOCATE (dumareaNA, STAT = DeAllocateStatus)
        IF (DeAllocateStatus /= 0) STOP "*** Deallocate failed ***"
      end if  
    end if

  DEALLOCATE (area_arr, STAT = DeAllocateStatus)
  IF (DeAllocateStatus /= 0) STOP "*** Deallocate failed ***"
  end if

end subroutine bourgouin

subroutine maskprecip3D(ptype,prect_mask,maskint)
  real, intent(inout) :: prect_mask(NLATS,NRECS)
  integer, intent(in) :: ptype(NLATS,NRECS)  
  integer, intent(in) :: maskint

  WHERE(ptype .EQ. maskint)
    prect_mask = prect_mask
  ELSEWHERE
    prect_mask = 0.0 
  END WHERE
end subroutine maskprecip3D








!$$$  Subprogram documentation block
!
! Subprogram: calwxt_bourg    Calculate precipitation type (Bourgouin)
!   Prgmmr: Baldwin      Org: np22        Date: 1999-07-06
!
! Abstract: This routine computes precipitation type
!    using a decision tree approach that uses the so-called
!    "energy method" of Bourgouin of AES (Canada) 1992
!
! Program history log:
!   1999-07-06  M Baldwin
!   1999-09-20  M Baldwin  make more consistent with bourgouin (1992)
!   2005-08-24  G Manikin  added to wrf post
!   2007-06-19  M Iredell  mersenne twister, best practices
!   2008-03-03  G Manikin  added checks to prevent stratospheric warming
!                           episodes from being seen as "warm" layers
!                           impacting precip type
!
! Usage:    call calwxt_bourg(im,jm,jsta_2l,jend_2u,jsta,jend,lm,lp1,   &
!    &                        iseed,g,pthresh,                          &
!    &                        t,q,pmid,pint,lmh,prec,zint,ptype)
!   Input argument list:
!     im       integer i dimension
!     jm       integer j dimension
!     jsta_2l  integer j dimension start point (including haloes)
!     jend_2u  integer j dimension end point (including haloes)
!     jsta     integer j dimension start point (excluding haloes)
!     jend     integer j dimension end point (excluding haloes)
!     lm       integer k dimension
!     lp1      integer k dimension plus 1
!     iseed    integer random number seed
!     g        real gravity (m/s**2)
!     pthresh  real precipitation threshold (m)
!     t        real(im,jsta_2l:jend_2u,lm) mid layer temp (K)
!     q        real(im,jsta_2l:jend_2u,lm) specific humidity (kg/kg)
!     pmid     real(im,jsta_2l:jend_2u,lm) mid layer pressure (Pa)
!     pint     real(im,jsta_2l:jend_2u,lp1) interface pressure (Pa)
!     lmh      real(im,jsta_2l:jend_2u) max number of layers
!     prec     real(im,jsta_2l:jend_2u) precipitation (m)
!     zint     real(im,jsta_2l:jend_2u,lp1) interface height (m)
!   Output argument list:
!     ptype    real(im,jm) instantaneous weather type ()
!              acts like a 4 bit binary
!                1111 = rain/freezing rain/ice pellets/snow
!                where the one's digit is for snow
!                      the two's digit is for ice pellets
!                      the four's digit is for freezing rain
!                  and the eight's digit is for rain
!              in other words...
!                ptype=1 snow
!                ptype=2 ice pellets/mix with ice pellets
!                ptype=4 freezing rain/mix with freezing rain
!                ptype=8 rain
!
! Modules used:
!   mersenne_twister pseudo-random number generator
!
! Subprograms called:
!   random_number    pseudo-random number generator
!
! Attributes:
!   Language: Fortran 90
!
! Remarks: vertical order of arrays must be layer   1 = top
!                                       and layer lmh = bottom
!
!$$$
      subroutine calwxt_bourg(rdseed,lm,lp1,t,pmid,pint,zint,ptype)

      implicit none
!
!    input:
      integer,intent(in):: lm,lp1,rdseed
      real,intent(in):: t(lm)
      real,intent(in):: pmid(lm)
      real,intent(in):: pint(lp1)
      real,intent(in):: zint(lp1)


!
!    output:
      integer,intent(out):: ptype
!
      integer :: ifrzl,iwrml,l,lhiwrm,lmhk
      real :: pintk1,areane,tlmhk,areape,pintk2,surfw,area1,dzkl,psfck
      real :: r1,r2
      real,parameter :: g=9.81

      integer, allocatable :: rndm_seed(:)
      integer :: rndm_seed_sz

!
!     initialize weather type array to zero (ie, off).
!     we do this since we want ptype to represent the
!     instantaneous weather type on return.



      ptype = 0

! Deal with random #'s

      call random_seed(size=rndm_seed_sz)
      allocate(rndm_seed(rndm_seed_sz))

      rndm_seed=rdseed
      call random_seed(put=rndm_seed)

      call random_number(r1)
      call random_number(r2)
      call random_number(r1)

      deallocate(rndm_seed)

! end

!
!      call random_number(rn,iseed)
!
!!$omp  parallel do
!!$omp& private(a,lmhk,tlmhk,iwrml,psfck,lhiwrm,pintk1,pintk2,area1,
!!$omp&         areape,dzkl,surfw,r1,r2)

      psfck=pint(lm+1)

!     find the depth of the warm layer based at the surface
!     this will be the cut off point between computing
!     the surface based warm air and the warm air aloft
!
!
!     lowest layer t
!
      tlmhk = t(lm)
      iwrml = lm + 1
      if (tlmhk.ge.273.15) then
        do l = lm, 2, -1
         if (t(l).ge.273.15.and.t(l-1).lt.273.15.and.           &
     &            iwrml.eq.lm+1) iwrml = l
          end do
      end if
!
!     now find the highest above freezing level
!
      lhiwrm = lm + 1
      do l = lm, 1, -1
! gsm  added 250 mb check to prevent stratospheric warming situations
!       from counting as warm layers aloft      
          if (t(l).ge.273.15 .and. pmid(l).gt.25000.) lhiwrm = l
      end do

!     energy variables
!     surfw is the positive energy between the ground
!     and the first sub-freezing layer above ground
!     areane is the negative energy between the ground
!     and the highest layer above ground
!     that is above freezing
!     areape is the positive energy "aloft"
!     which is the warm energy not based at the ground
!     (the total warm energy = surfw + areape)
!
!     pintk1 is the pressure at the bottom of the layer
!     pintk2 is the pressure at the top of the layer
!     dzkl is the thickness of the layer
!     ifrzl is a flag that tells us if we have hit
!     a below freezing layer
!
      pintk1 = psfck
      ifrzl = 0
      areane = 0.0
      areape = 0.0
      surfw = 0.0                                         

      do l = lm, 1, -1
          if (ifrzl.eq.0.and.t(l).le.273.15) ifrzl = 1
          pintk2=pint(l)
          dzkl=zint(l)-zint(l+1)
          area1 = log(t(l)/273.15) * g * dzkl
          if (t(l).ge.273.15.and. pmid(l).gt.25000.) then
              if (l.lt.iwrml) areape = areape + area1
              if (l.ge.iwrml) surfw = surfw + area1
          else
              if (l.gt.lhiwrm) areane = areane + abs(area1)
          end if
          pintk1 = pintk2
      end do
      
!
!     decision tree time
!
      if (areape.lt.2.0) then
!         very little or no positive energy aloft, check for
!         positive energy just above the surface to determine rain vs. snow
          if (surfw.lt.5.6) then
!             not enough positive energy just above the surface
!             snow = 0
              ptype = 0
          else if (surfw.gt.13.2) then
!             enough positive energy just above the surface
!             rain = 2
              ptype = 2
          else
!             transition zone, assume equally likely rain/snow
!             picking a random number, if <=0.5 snow
              !r1 = rn(1)
              if (r1.le.0.5) then
!                 snow = 0
                  ptype = 0
              else
!                 rain = 2
                  ptype = 2
              end if
          end if
!
      else
!         some positive energy aloft, check for enough negative energy
!         to freeze and make ice pellets to determine ip vs. zr
          if (areane.gt.66.0+0.66*areape) then
!             enough negative area to make ip,
!             now need to check if there is enough positive energy
!             just above the surface to melt ip to make rain
              if (surfw.lt.5.6) then
!                 not enough energy at the surface to melt ip
!                 ice pellets = 3
                  ptype = 3
              else if (surfw.gt.13.2) then
!                 enough energy at the surface to melt ip
!                 rain = 2
                  ptype = 2
              else
!                 transition zone, assume equally likely ip/rain
!                 picking a random number, if <=0.5 ip
                  !r1 = rn(1)
                  if (r1.le.0.5) then
!                     ice pellets = 3
                      ptype = 3
                  else
!                     rain = 2
                      ptype = 2
                  end if
              end if
          else if (areane.lt.46.0+0.66*areape) then
!             not enough negative energy to refreeze, check surface temp
!             to determine rain vs. zr
              if (tlmhk.lt.273.15) then
!                 freezing rain = 4
                  ptype = 4
              else
!                 rain = 2
                  ptype = 2
              end if
          else
!             transition zone, assume equally likely ip/zr
!             picking a random number, if <=0.5 ip
              !r1 = rn(1)
              if (r1.le.0.5) then
!                 still need to check positive energy
!                 just above the surface to melt ip vs. rain
                  if (surfw.lt.5.6) then
!                     ice pellets = 3
                      ptype = 3
                  else if (surfw.gt.13.2) then
!                     rain = 2
                      ptype = 2
                  else
!                     transition zone, assume equally likely ip/rain
!                     picking a random number, if <=0.5 ip
                      !r2 = rn(2)
                      if (r2.le.0.5) then
!                         ice pellets = 3
                          ptype = 3
                      else
!                         rain = 2
                          ptype = 2
                      end if
                  end if
              else
!                 not enough negative energy to refreeze, check surface temp
!                 to determine rain vs. zr
                  if (tlmhk.lt.273.15) then
!                     freezing rain = 4
                      ptype = 4
                  else
!                     rain = 2
                      ptype = 2
                  end if
              end if
          end if
      end if
!      end do
!      end do
      return
end subroutine calwxt_bourg




subroutine init_random_seed

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)

end subroutine init_random_seed



end program pres_temp_4D_rd

