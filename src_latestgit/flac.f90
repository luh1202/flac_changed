!-*- F90 -*-

! --------- Flac ------------------------- 

subroutine flac

use arrays
include 'precision.inc' 
include 'params.inc'
include 'arrays.inc'

! Update Thermal State
! Skip the therm calculations if itherm = 3
call fl_injectheat
call fl_therm

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

! Changing marker phases
! XXX: change_phase is slow, don't call it every loop
!Hao commented following line to skip phase change 5.15.2018
if( mod(nloop, 10).eq.0 ) call change_phase_dike
!if( mod(nloop, 100).eq.0 ) call write_stress
! Update stresses by constitutive law (and mix isotropic stresses)
call fl_rheol
!print*,'hey there'
!call user_lu
!end if
!function write_stress!(iz,ix)
!use arrays
!include 'precision.inc'
!include 'params.inc'
!include 'arrays.inc'
!s11 = 0.25 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
!if ((iz==2).and.(ix==44)) then
!open (unit = 1, file = "s111.txt")
!write (1,*) "Here are the s111 ", s11
!close (1)
!!end if
!!s22 = 0.25 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
!!if ((iz=2).and.(ix=44)) then
!open (unit = 1, file = "s122.txt")
!write (1,*) "Here are the s122 ", s22
!close (1)
!!end if
!!s33 = 0.25 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
!!if ((iz=2).and.(ix=44)) then
!open (unit = 1, file = "s133.txt")
!write (1,*) "Here are the s133 ", s33
!close (1)
!end if

! update stress boundary conditions
if (ynstressbc.eq.1.) call bc_update

! Calculations in a node: forces, balance, velocities, new coordinates
call fl_node
!call functions
!if( mod(nloop, 100).eq.0 ) call write_stress

! New coordinates
call fl_move

! Adjust real masses due to temperature
if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses

! Adjust inertial masses or time step due to deformations
if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass

! Adjust time Step 
call dt_adjust

500 continue


return
end subroutine flac
