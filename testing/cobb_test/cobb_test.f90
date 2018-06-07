program main

!*****************************************************************************
! This is a program to test the Cobb subroutine written for CAM
! Colin Zarzycki 7/9/15
!*****************************************************************************

  use mod_cobb

  implicit none

  real,dimension(3) :: temp  = (/248.0,258.0,268.0/)
  real,dimension(3) :: omega = (/5.0,8.0,12.0/)
  real,dimension(3) :: rh = (/95.0,95.0,95.0/)
  real,dimension(3) :: pres = (/85000.,90000.,95000./)
  real,dimension(4) :: z = (/2000.,1500.,1000.,500./)
  real,dimension(3) :: zm = (/2000.,1500.,1000./)

  real :: snowratio

  real :: r1,r2

  snowratio = 0.0

  call init_random_seed()

  call random_number(r2)
  call random_number(r1)
  call random_number(r2)

  print *,r1
  print *,r2


  call cobb(omega,temp,rh,z,pres,snowratio)

  print *,snowratio

  stop
end
