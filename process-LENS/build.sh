#!/bin/bash

ifort -fPIC -shared-intel -mcmodel=large process_snow.f90 omcalc_ccm.f rvdv.f phybrid_ccm.f int2p_dp.f relhum_dp.f -o 1out

