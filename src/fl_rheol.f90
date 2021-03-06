
!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  

subroutine fl_rheol
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension depl(4)
dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
logical rh_sel
pi = 3.1415926
!if( mod(nloop,10).eq.0 .OR. ireset.eq.1 ) then
!    rh_sel = .true.
!else
!    rh_sel = .false.
!endif
rh_sel = .true.

!XXX: irh==11, or irh>=11?
irh=irheol(mphase)
if(irh.eq.11) call init_visc

!if(iynts.eq.1) call init_temp

! Initial stress boundary condition
! Accretional Stresses
if (ny_inject.gt.0) then
  sarc1 = 0.
  sarc2 = 0. 
  if (ny_inject.eq.1) iinj = 1
  if (ny_inject.eq.2) iinj = (nx-1)/2

  dxinj = 0.
  do jinj = 1,nelem_inject
    dxinj=dxinj+cord(jinj,iinj+1,1)-cord(jinj,iinj,1)
  enddo
  dxinj = dxinj/nelem_inject
endif

irh_mark = 0

! max. deviatoric strain and area change of current time step
curr_devmax = devmax
curr_dvmax = dvmax

!$OMP Parallel Private(i,j,k,iph,irh,bulkm,rmu,coh,phi,psi, &
!$OMP                  stherm,hardn,vis, &
!$OMP                  de11,de22,de12,de33,dv, &
!$OMP                  s11p,s22p,s12p,s33p, &
!$OMP                  s11v,s22v,s12v,s33v, &
!$OMP                  depl,ipls,diss, &
!$OMP                  sII_plas,sII_visc, &
!$OMP                  quad_area,s0a,s0b,s0) &
!$OMP firstprivate(irh_mark)
!$OMP do schedule(guided) reduction(max: curr_devmax, curr_dvmax)
do 3 i = 1,nx-1 !nx-1 = 120 element number on x direction
    do 3 j = 1,nz-1 !nz-1 = 40 element number on z direction
        ! iphase (j,i) is number of a phase NOT a rheology
        !print *, 'nx-1 =', nx-1
        !print *, 'nz-1 =', nz-1
        iph = iphase(j,i)
        irh = irheol(iph)
        zcord_ave = 0.25 * (cord(j,i,2) + cord(j+1,i,2) + cord(j,i+1,2) + cord(j+1,i+1,2))
        temp_ave = 0.25 * (temp(j,i) + temp(j+1,i) + temp(j,i+1) + temp(j+1,i+1))
!        if (it.eq.1) then
!         if (temp_ave.le.600) then
!             rate_inject = rate_inject_brittle
!         elseif (temp_ave.gt.600) then
!             rate_inject = rate_inject_ductile
!         endif
!        endif

        !if (it.eq.2) then
         if (temp_ave.le.600) then
             rate_inject = rate_inject_brittle
         elseif ((temp_ave.gt.600).and.(time .lt. time_max*0.3)) then
             rate_inject = rate_inject_ductile

         elseif ((temp_ave.gt.600).and.(time .ge. time_max*0.3).and.(time .lt. time_max*0.6)) then
             rate_inject = rate_inject_ductile_e

         elseif ((temp_ave.gt.600).and.(time .ge. time_max*0.6).and.(time .lt. time_max)) then
             rate_inject = rate_inject_ductile_s

             !rate_inject = (fa * SIN(2 * pi/fb * time) + fc) * (fsr)  !sine function
             !print*,'fsr=',fsr
             !print*, 'rate =', rate_inject
             !print*,'it=',it
         endif
        !endif
!print*, 'it =', it
    if(ny_inject.gt.0.and.j.le.nelem_inject) then
!!!!!!!!!!!!!!!!!!!!!!!!!change in rigidity!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!too high!!!!!!!!!!!!!!
!        if (it.eq.1) then
!        if (temp_ave.lt.600. .and. (zcord_top-zcord_ave) .le. 7e3 .and.(fa * SIN(2 * pi/fb * time) + fc).ge.0.6) then
!
!           poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
!           young = rm(iph)*2.*(1.+poiss)
!           sarc1 = -young/(1.-poiss*poiss)*rate_inject/dxinj*dt
!           sarc2 = sarc1*poiss/(1.-poiss)
!        !!!!!!!!!!!too low!!!!!!!!!!!!!!!
!!        elseif (temp_ave.lt.600. .and. (zcord_top-zcord_ave) .le. 7e3 .and.(fa * SIN(2 * pi/fb * time) + fc).le.0.25) then
!!           poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
!!           young = 10. * rm(iph)*2.*(1.+poiss)
!!           sarc1 = -young/(1.-poiss*poiss)*rate_inject/dxinj*dt
!!           sarc2 = sarc1*poiss/(1.-poiss)
!          else

!           poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
!           young = rm(iph)*2.*(1.+poiss)
!           sarc1 = -young/(1.-poiss*poiss)*rate_inject/dxinj*dt
!           sarc2 = sarc1*poiss/(1.-poiss)
!          endif
!        endif
!!!!!!!!!!!!!!!!!!!!!!!!change in rigidity ends!!!!!!!!!!!!!!!!!!
        !if (it.eq.2) then
         poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
         young = rm(iph)*2.*(1.+poiss)
         sarc1 = -young/(1.-poiss*poiss)*rate_inject/dxinj*dt
         sarc2 = sarc1*poiss/(1.-poiss)
        endif
     !endif
        bulkm = rl(iph) + 2.*rm(iph)/3.
        rmu   = rm(iph)

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))

        ! Preparation of plastic properties
        !if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)
if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)

        ! Re-evaluate viscosity
        if (irh.eq.7 .or. irh.eq.12) then
!if (irh.eq.7 .or. irh.eq.12) then
            if( mod(nloop,ifreq_visc).eq.0 .OR. ireset.eq.1 ) visn(j,i) = Eff_visc(j,i)
        !            if (ny_inject.gt.0.and.i.eq.iinj) visn(j,i) = v_min
        endif
        vis = visn(j,i)
        !print *,'6'
        ! Cycle by triangles
        do k = 1,4
            ! Incremental strains
            de11 = strainr(1,k,j,i)*dt
            de22 = strainr(2,k,j,i)*dt
            de12 = strainr(3,k,j,i)*dt
            de33 = 0.
            dv = dvol(j,i,k)
            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
            if(ny_inject.gt.0.and.j.le.nelem_inject) then
                !XXX: iinj is un-init'd if ny_inject is not 1 or 2.
                if(i.eq.iinj) then
                    s11p(k) = stress0(j,i,1,k) + sarc1
                    s22p(k) = stress0(j,i,2,k) + sarc2
                    !!            irh = 1
                endif
            endif
            !print *,'7'
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)
        !!            if(abs(sarc11).gt.0.) write(*,*) i,j,sarc11,sarc22
                if (irh.eq.1) then
!if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                irheol_fl(j,i) = 0  
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)
            !print *,'8'
                elseif (irh.eq.7) then
!elseif (irh.eq.7) then
                ! viscous
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,&
                    ndim,dt,curr_devmax,curr_dvmax)
                irheol_fl(j,i) = -1  
                stress0(j,i,1,k) = s11v(k)
                stress0(j,i,2,k) = s22v(k)
                stress0(j,i,3,k) = s12v(k)
                stress0(j,i,4,k) = s33v(k)
                elseif (irh.eq.6) then
            !elseif (irh.eq.6) then
                ! plastic
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                    ten_off,ndim,irh_mark)
                irheol_fl(j,i) = 1
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)
            !print *,'9'
             elseif (irh.ge.11) then
            !elseif (irh.ge.11) then
                ! Mixed rheology (Maxwell or plastic)
                if( rh_sel ) then
                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                        s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                        ten_off,ndim,irh_mark)
                    call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                        de11,de22,de33,de12,dv,&
                        ndim,dt,curr_devmax,curr_dvmax)
                else ! use previously defined rheology
                    if( irheol_fl(j,i) .eq. 1 ) then
                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                            s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                            ten_off,ndim,irh_mark)
                        stress0(j,i,1,k) = s11p(k)
                        stress0(j,i,2,k) = s22p(k)
                        stress0(j,i,3,k) = s12p(k)
                        stress0(j,i,4,k) = s33p(k)
                    else  ! irheol_fl(j,i) = -1
                        call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                            de11,de22,de33,de12,dv,&
                            ndim,dt,curr_devmax,curr_dvmax)
                        stress0(j,i,1,k) = s11v(k)
                        stress0(j,i,2,k) = s22v(k)
                        stress0(j,i,3,k) = s12v(k)
                        stress0(j,i,4,k) = s33v(k)
                    endif
                endif
            endif
        enddo
        !print *,'10'
        if( irh.ge.11 .AND. rh_sel ) then
            ! deside - elasto-plastic or viscous deformation
            sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                    + 4*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2

            sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                    + 4*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

            if (sII_plas .lt. sII_visc) then
                do k = 1, 4
                    stress0(j,i,1,k) = s11p(k)
                    stress0(j,i,2,k) = s22p(k)
                    stress0(j,i,3,k) = s12p(k)
                    stress0(j,i,4,k) = s33p(k)
                end do
                irheol_fl (j,i) = 1
            else 
                do k = 1, 4
                    stress0(j,i,1,k) = s11v(k)
                    stress0(j,i,2,k) = s22v(k)
                    stress0(j,i,3,k) = s12v(k)
                    stress0(j,i,4,k) = s33v(k)
                end do
                irheol_fl (j,i) = -1
            endif
        endif
        !print *,'11'

        ! Averaging of isotropic stresses for pair of elements
        if (mix_stress .eq. 1 ) then

            ! For A and B couple:
            ! area(n,it) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1./(area(j,i,1)+area(j,i,2))
            s0a=0.5*(stress0(j,i,1,1)+stress0(j,i,2,1))
            s0b=0.5*(stress0(j,i,1,2)+stress0(j,i,2,2))
            s0=(s0a*area(j,i,2)+s0b*area(j,i,1))*quad_area
            stress0(j,i,1,1) = stress0(j,i,1,1) - s0a + s0
            stress0(j,i,2,1) = stress0(j,i,2,1) - s0a + s0
            stress0(j,i,1,2) = stress0(j,i,1,2) - s0b + s0
            stress0(j,i,2,2) = stress0(j,i,2,2) - s0b + s0

            ! For C and D couple:
            quad_area = 1./(area(j,i,3)+area(j,i,4))
            s0a=0.5*(stress0(j,i,1,3)+stress0(j,i,2,3))
            s0b=0.5*(stress0(j,i,1,4)+stress0(j,i,2,4))
            s0=(s0a*area(j,i,4)+s0b*area(j,i,3))*quad_area
            stress0(j,i,1,3) = stress0(j,i,1,3) - s0a + s0
            stress0(j,i,2,3) = stress0(j,i,2,3) - s0a + s0
            stress0(j,i,1,4) = stress0(j,i,1,4) - s0b + s0
            stress0(j,i,2,4) = stress0(j,i,2,4) - s0b + s0
            !print *, 'stress0(2,44,1,3)=', stress0(2,44,1,3)
!            sxx = 0.25 * (stress0(2,44,1,1)+stress0(2,44,1,2)+stress0(2,44,1,3)+stress0(2,44,1,4))
!            syy = 0.25 * (stress0(2,44,2,1)+stress0(2,44,2,2)+stress0(2,44,2,3)+stress0(2,44,2,4))
!            szz = 0.25 * (stress0(2,44,4,1)+stress0(2,44,4,2)+stress0(2,44,4,3)+stress0(2,44,4,4))
!                !if((jj==2).and.(ii==44)) then
!                    open (unit = 1, file = "sxx.txt")
!                        write (1,*) "Here are the sxx ", sxx
!                    !close (1)
!
!                    open (unit = 2, file = "syy.txt")
!                        write (2,*) "Here are the syy ", syy
!                    !close (1)
!
!                    open (unit = 3, file = "szz.txt")
!                        write (3,*) "Here are the szz ", szz
                    !close (1)
                !end if
!                end do
!            end do

        endif
        !print *,'12'
        if (irh.eq.6 .or. irh.ge.11) then
            !  ACCUMULATED PLASTIC STRAIN
            ! Average the strain for pair of the triangles
            ! Note that area (n,it) is inverse of double area !!!!!
            dps = 0.5*( depl(1)*area(j,i,2)+depl(2)*area(j,i,1) ) / (area(j,i,1)+area(j,i,2)) &
                + 0.5*( depl(3)*area(j,i,4)+depl(4)*area(j,i,3) ) / (area(j,i,3)+area(j,i,4))
            aps(j,i) = aps(j,i) + dps 
            if( aps(j,i) .lt. 0. ) aps(j,i) = 0.

            !	write(*,*) depl(1),depl(2),depl(3),depl(4),area(j,i,1),area(j,i,2),area(j,i,3),area(j,i,4)
            !print *,'13'
            ! LINEAR HEALING OF THE PLASTIC STRAIN
            !if (tau_heal .ne. 0. .and. dps .le. 0.) then
            ! New healing by Eunseo Choi, Jul 2018
            if (tau_heal .ne. 0. .and. aps(j,i) .gt. 0.) then
              sr11 = 0.25 * (strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
              sr22 = 0.25 * (strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
              sr12 = 0.25 * (strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))
              srJ2 = 0.5 * sqrt((sr11-sr22)**2 + 4*sr12*sr12)
              daps = -aps(j,i)/tau_heal
              if (srJ2 .ge. 1.e-13) then
               daps = daps + srJ2
              endif
              aps(j,i) = aps(j,i) + dt * daps
              if( aps(j,i) .lt. 0. ) aps(j,i) = 0.
            endif
            if (ny_inject.gt.0 .and. i.eq.iinj) aps (j,i) = 0.
            !if (ny_inject.gt.0.and. (i.le.(iinj+1))) aps (j,i) = 0.
        end if


!!!!!!!!!!!!debug!!!!!!!!!!!!!
sxx5 = 0.25 * (stress0(5,26,1,1)+stress0(5,26,1,2)+stress0(5,26,1,3)+stress0(5,26,1,4))
syy5 = 0.25 * (stress0(5,26,2,1)+stress0(5,26,2,2)+stress0(5,26,2,3)+stress0(5,26,2,4))
szz5 = 0.25 * (stress0(5,26,4,1)+stress0(5,26,4,2)+stress0(5,26,4,3)+stress0(5,26,4,4))
sxx4 = 0.25 * (stress0(4,26,1,1)+stress0(4,26,1,2)+stress0(4,26,1,3)+stress0(4,26,1,4))
syy4 = 0.25 * (stress0(4,26,2,1)+stress0(4,26,2,2)+stress0(4,26,2,3)+stress0(4,26,2,4))
szz4 = 0.25 * (stress0(4,26,4,1)+stress0(4,26,4,2)+stress0(4,26,4,3)+stress0(4,26,4,4))
sxx6 = 0.25 * (stress0(6,26,1,1)+stress0(6,26,1,2)+stress0(6,26,1,3)+stress0(6,26,1,4))
syy6 = 0.25 * (stress0(6,26,2,1)+stress0(6,26,2,2)+stress0(6,26,2,3)+stress0(6,26,2,4))
szz6 = 0.25 * (stress0(6,26,4,1)+stress0(6,26,4,2)+stress0(6,26,4,3)+stress0(6,26,4,4))

sxx9 = 0.25 * (stress0(9,62,1,1)+stress0(9,62,1,2)+stress0(9,62,1,3)+stress0(9,62,1,4))
syy9 = 0.25 * (stress0(9,62,2,1)+stress0(9,62,2,2)+stress0(9,62,2,3)+stress0(9,62,2,4))
szz9 = 0.25 * (stress0(9,62,4,1)+stress0(9,62,4,2)+stress0(9,62,4,3)+stress0(9,62,4,4))
sxx8 = 0.25 * (stress0(8,62,1,1)+stress0(8,62,1,2)+stress0(8,62,1,3)+stress0(8,62,1,4))
syy8 = 0.25 * (stress0(8,62,2,1)+stress0(8,62,2,2)+stress0(8,62,2,3)+stress0(8,62,2,4))
szz8 = 0.25 * (stress0(8,62,4,1)+stress0(8,62,4,2)+stress0(8,62,4,3)+stress0(8,62,4,4))
sxx10 = 0.25 * (stress0(10,62,1,1)+stress0(10,62,1,2)+stress0(10,62,1,3)+stress0(10,62,1,4))
syy10 = 0.25 * (stress0(10,62,2,1)+stress0(10,62,2,2)+stress0(10,62,2,3)+stress0(10,62,2,4))
szz10 = 0.25 * (stress0(10,62,4,1)+stress0(10,62,4,2)+stress0(10,62,4,3)+stress0(10,62,4,4))
sec_year = 3.1558e+7

if( mod(nloop, 50000).eq.0 ) then
open (unit = 1, file = "sxx5.txt",position='append')
open (unit = 2, file = "syy5.txt",position='append')
open (unit = 3, file = "szz5.txt",position='append')
write (1,*) sxx5
write (2,*) syy5
write (3,*) szz5

open (unit = 4, file = "sxx4.txt",position='append')
open (unit = 5, file = "syy4.txt",position='append')
open (unit = 6, file = "szz4.txt",position='append')
write (4,*) sxx4
write (5,*) syy4
write (6,*) szz4

open (unit = 7, file = "sxx6.txt",position='append')
open (unit = 8, file = "syy6.txt",position='append')
open (unit = 9, file = "szz6.txt",position='append')
write (7,*) sxx6
write (8,*) syy6
write (9,*) szz6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open (unit = 10, file = "sxx9.txt",position='append')
open (unit = 11, file = "syy9.txt",position='append')
open (unit = 12, file = "szz9.txt",position='append')
write (10,*) sxx9
write (11,*) syy9
write (12,*) szz9

open (unit = 13, file = "sxx8.txt",position='append')
open (unit = 14, file = "syy8.txt",position='append')
open (unit = 15, file = "szz8.txt",position='append')
write (13,*) sxx8
write (14,*) syy8
write (15,*) szz8

open (unit = 16, file = "sxx10.txt",position='append')
open (unit = 17, file = "syy10.txt",position='append')
open (unit = 18, file = "szz10.txt",position='append')
write (16,*) sxx10
write (17,*) syy10
write (18,*) szz10

open (unit = 19, file = "timeinyr.txt",position='append')
write (19,*) (time/sec_year/1.e6)

end if
!!!!!!!!!!!debug ends!!!!!!!!!!!!!!


        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain(j,i,2) = strain(j,i,2) + 0.25*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain(j,i,3) = strain(j,i,3) + 0.25*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))
!!!!!!!!!!!!!!!!!!!!!!!debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        strain_1st = 0.5 * ( strain(5,26,1) + strain(5,26,2) )
!        strain_2nd = 0.5 * sqrt((strain(5,26,1)-strain(5,26,2))**2 + 4*strain(5,26,3)**2)
!        if( mod(nloop, 50000).eq.0 ) then
!        open (unit = 7, file = "strain_1st.txt")
!        open (unit = 8, file = "strain_2nd.txt")
!        write(7,*) strain_1st
!        write (8,*) strain_2nd
!        end if
!!!!!!!!!!!!!!!!!!!!!!debug ends!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print *,'strain1 =', strain(2,44,1)
        !print *,'strain2 =', strain(2,44,2)
        !print *,'strain3 =', strain(2,44,3)
!do  i = 1,nx-1 !nx-1 = 120 element number on x direction
!do  j = 1,nz-1

!end do
!end do

!!!!!!!!!!!!!!!!!!!!debug starts!!!!!!!!!!!!!!!!!!!
!do  ii = 1,nx-1 !nx-1 = 120 element number on x direction
!do  jj = 1,nz-1
!iph = iphase(jj,ii)
!irh = irheol(iph)
!do k = 1,4
!if((jj==5).and.(ii==26)) then
!call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
!ten_off,ndim,irh_mark)
!open (unit = 1, file = "s11.txt")
!write (1,*) s11
!
!open (unit = 2, file = "s22.txt")
!write (2,*) s22
!
!open (unit = 3, file = "anphi.txt")
!write (3,*) anphi
!
!open (unit = 4, file = "amc.txt")
!write (4,*) amc
!
!end if
!enddo
!enddo
!enddo
!!!!!!!!!!!!!!!!!!!!debug ends!!!!!!!!!!!!!!!!!!!!
3 continue
!$OMP end do
!$OMP end parallel


devmax = max(devmax, curr_devmax)
dvmax = max(dvmax, curr_dvmax)

return
end
