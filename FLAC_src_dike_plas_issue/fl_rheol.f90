!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  

subroutine fl_rheol
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
dimension njTinj(nz)
dimension depl(4)
dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
logical rh_sel
pi = 3.1415926

rh_sel = .true.

!XXX: irh==11, or irh>=11?
irh=irheol(mphase)
if(irh.ge.11) call init_visc
if(iynts.eq.1) call init_temp

! Initial stress boundary condition
! Accretional Stresses
if (ny_inject.gt.0) then
sarc1 = 0.
sarc2 = 0.
if (ny_inject.eq.1) iinj = 1
if (ny_inject.eq.2) iinj = nx/2
dxinj = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!original sarc1&2 determaination part!!!!!!!!!!!!!!!!!!!!!!
do jinj = 1,nelem_inject
iph=iphase(jinj,iinj)
dxinj=dxinj+cord(jinj,iinj+1,1)-cord(jinj,iinj,1)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!
dxinj = dxinj/nelem_inject
!Constants Elastic:
poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
young = rm(iph)*2.*(1.+poiss)
!Additional Stress:
sarc1 = -young/(1.-poiss*poiss)*rate_inject/dxinj*dt
sarc2 = sarc1*poiss/(1.-poiss)
endif
!!!!!!!!!!!!!!!!new!!!!!!!!!!!!!!!!!!
inj_count = inj_count + 1
!!!!!!!!!!!!!!!ends!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!new!!!!!!!!!!!!!!!
ninjtop=jinj1
ninjbot=jinj2
ninj=iinj2-iinj1+1
nelem_inject=jinj2-jinj1+1
!!!!!!!!!!!!!!!!!new ends!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!new!!!!!!!!!!!!!!!!!!!
   !rogh = 0. !Reset for pressure calculation [M. Behn April 4, 2006]
! Calculate ninjbot to be the minimum of ninjbot and the depth to the maximum cutoff temp (Tinj)
!   jcnt=1
!   njTinj(1) = ninjbot
!    do j = 1,nz-1
!        !dcord = cord(1,i,2) - cord(j+1,i,2)
!        dcord = cord(j+1,i,2) - cord(1,i,2)
!         if(dcord.gt.Tcinj) then
!           njTinj(jcnt) = j
!           jcnt = jcnt+1
!         endif
!    end do
!    ninjbot = min(ninjbot,njTinj(1))
!!!!!!!!!!!!!!!!!!new ends!!!!!!!!!!!!!!!!!!!!!
    do 3 j = 1,nz-1 !nz-1 = 40 element number on z direction
        ! iphase (j,i) is number of a phase NOT a rheology
        iph = iphase(j,i)
        irh = irheol(iph)
        zcord_ave = 0.25 * (cord(j,i,2) + cord(j+1,i,2) + cord(j,i+1,2) + cord(j+1,i+1,2))
        temp_ave = 0.25 * (temp(j,i) + temp(j+1,i) + temp(j,i+1) + temp(j+1,i+1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!new starts!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        densT = den(iph) * ( 1 - alfa(iph)*temp_ave )
!        dh1 = cord (j,i  ,2) - cord (j+1,i  ,2)
!        dh2 = cord (j,i+1,2) - cord (j+1,i+1,2)
!        dh  = 0.5 * (dh1+dh2)
!        dPT = densT * g * dh
!        dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )
!        press = rogh + 0.5*dP
!        rogh = rogh + dP
!!!!!!!!!!!!!!!!!!!!!!!!!!new ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if (it.eq.1) then
         if (temp_ave.le.600) then
             rate_inject = rate_inject_brittle
         elseif (temp_ave.gt.600) then
             rate_inject = rate_inject_ductile
         endif
        !endif

        !if (it.eq.2) then
!         if (temp_ave.le.600) then
!             rate_inject = rate_inject_brittle
!         elseif ((temp_ave.gt.600).and.(time .lt. time_max*0.3)) then
!             rate_inject = rate_inject_ductile
!
!         elseif ((temp_ave.gt.600).and.(time .ge. time_max*0.3).and.(time .lt. time_max*0.6)) then
!             rate_inject = rate_inject_ductile_e
!
!         elseif ((temp_ave.gt.600).and.(time .ge. time_max*0.6).and.(time .lt. time_max)) then
!             rate_inject = rate_inject_ductile_s
!         endif
    !nelem_inject = int(dike_depth/(ABS(rzbo/(nz-1))))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!original starts!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ny_inject.gt.0 .and. j.le.nelem_inject) then
         poiss = 0.5 * rl(iph) / (rl(iph)+rm(iph))
         young = rm(iph)*2.*(1. + poiss)
         fdum = ratfac*rate_inject*dt/(dxinj*dble(ninj))
         sarc1 = -young / (1. - (poiss*poiss)) * fdum
         sarc2 = sarc1*poiss/(1. - poiss)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!original ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!new ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bulkm = rl(iph) + 2.*rm(iph)/3.
rmu   = rm(iph)
coh   = coha(iph)
phi   = phimean(iph)
psi   = psia(iph)
poiss = 0.5 * rl(iph) / (rl(iph)+rm(iph))
young = rm(iph)*2.*(1. + poiss)

if ((ny_inject.le.2)) then
if ((i.ge.iinj1).and.(i.le.iinj2).and.(j.le.ninjbot).and.(j.ge.ninjtop)) then
dxinj = 0.5*(cord(j,i+1,1)-cord(j,i,1)+cord(j+1,i+1,1)-cord(j+1,i,1))
fdum = ratfac*rate_inject*dt/(dxinj*dble(ninj))
sarc1v = -young/(1.-poiss*poiss)*fdum
sarc2v = poiss*sarc1v
sarc3v = 0.

sarc1p = -(1.-poiss)*young/(1.+poiss)/(1.-2.*poiss)*fdum     ! Plane Strain Dikes - M. Behn 04/04/05
sarc2p = sarc1p*poiss/(1.-poiss)
sarc3p = poiss * (sarc1p + sarc2p)
endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!new ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))

        ! Preparation of plastic properties
        if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)

        ! Re-evaluate viscosity
        if (irh.eq.7 .or. irh.eq.12) then
            if( mod(nloop,ifreq_visc).eq.0 .OR. ireset.eq.1 ) visn(j,i) = Eff_visc(j,i)
        !            if (ny_inject.gt.0.and.i.eq.iinj) visn(j,i) = v_min
        endif
        !endif
        vis = visn(j,i)
        ! Cycle by triangles
        do  k = 1,4
            ! Incremental strains
            de11 = strainr(1,k,j,i)*dt
            de22 = strainr(2,k,j,i)*dt
            de12 = strainr(3,k,j,i)*dt
            de33 = 0.
            dv = dvol(j,i,k)
            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
!            if(ny_inject.gt.0.and.j.le.nelem_inject) then
!                !XXX: iinj is un-init'd if ny_inject is not 1 or 2.
!                if(i.eq.iinj) then
!!!!!!!!!!!!!!!!!!!!!!!!!!original starts!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    s11p(k) = stress0(j,i,1,k) + sarc1
!!                    s22p(k) = stress0(j,i,2,k) + sarc2
!!!!!!!!!!!!!!!!!!!!!!!!!original ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!new starts!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    s11p(k) = stress0(j,i,1,k)  + stherm
!                    s22p(k) = stress0(j,i,2,k)  + stherm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!new ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    !!            irh = 1
!                endif
!            endif
            !print *,'7'
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)
        !!            if(abs(sarc11).gt.0.) write(*,*) i,j,sarc11,sarc22
                if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                irheol_fl(j,i) = 0  
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)
            !print *,'8'
                elseif (irh.eq.7) then
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
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                    s11p(k),s22p(k),s33p(k),s12p(k),&
                    de11,de22,de33,de12,ten_off,ndim,irh_mark)

                irheol_fl(j,i) = 1
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)
            !print *,'9'
             elseif (irh.ge.11) then
            !elseif (irh.ge.11) then
                ! Mixed rheology (Maxwell or plastic)
                if( rh_sel ) then !always the case
                depl(k) = 0.
                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                        s11p(k),s22p(k),s33p(k),s12p(k),&
                        de11,de22,de33,de12,ten_off,ndim,irh_mark)

                    call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                        de11,de22,de33,de12,dv,&
                        ndim,dt,curr_devmax,curr_dvmax)
                else ! use previously defined rheology
                    if( irheol_fl(j,i) .eq. 1 ) then
                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                            s11p(k),s22p(k),s33p(k),s12p(k),&
                            de11,de22,de33,de12,ten_off,ndim,irh_mark)

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
            !print *, 'irheol_fl(j,i) =',irheol_fl(j,i) 1 and -1 show up mixed
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
!!!!!!!!!!!!!!!!!!!!!starts!!!!!!!!!!!!!!!!!

if((ny_inject.eq.2).and.(inj_count.gt.ext_period)) then
if((i.ge.iinj1 ).and.(i.le.iinj2).and.(j.le.ninjbot).and.(j.ge.ninjtop)) then
!if(i.eq.iinj-2) then

ntap = 0
if(ntap.eq.0) then
dtpr=1.
elseif(ntap.eq.1) then
dtpr = (j + 0.5 - ninjtop) * (2./nelem_inject)
elseif(ntap.eq.2) then
dtpr = sqrt( (nelem_inject/2.)**2 - ((j-0.5) - ((ninjtop+ninjbot-1.)/2.))**2 ) / (nelem_inject/2.)
endif

stress0(j,i,1,k) = s11p(k) +sarc1p*dtpr
stress0(j,i,2,k) = s22p(k) +sarc2p*dtpr
stress0(j,i,4,k) = s33p(k) +sarc3p*dtpr  ! Reset s33 for accretion M. Behn 04/19/05
endif

endif
!!!!!!!!!!!!!!!!!!!!!ends!!!!!!!!!!!!!!!!!!!
                end do
                irheol_fl (j,i) = 1
            else 
                do k = 1, 4
                    stress0(j,i,1,k) = s11v(k)
                    stress0(j,i,2,k) = s22v(k)
                    stress0(j,i,3,k) = s12v(k)
                    stress0(j,i,4,k) = s33v(k)
!!!!!!!!!!!!!!!!!!!!!starts!!!!!!!!!!!!!!!!!
if((ny_inject.eq.2).and.(inj_count.gt.ext_period)) then
if((i.ge.iinj1 ).and.(i.le.iinj2).and.(j.le.ninjbot).and.(j.ge.ninjtop)) then
!if(i.eq.iinj-2) then

ntap = 0
if(ntap.eq.0) then
dtpr=1.
elseif(ntap.eq.1) then
dtpr = (j + 0.5 - ninjtop) * (2./nelem_inject)
elseif(ntap.eq.2) then
dtpr = sqrt( (nelem_inject/2.)**2 - ((j-0.5) - ((ninjtop+ninjbot-1.)/2.))**2 ) / (nelem_inject/2.)
endif

stress0(j,i,1,k) = s11v(k) +sarc1v*dtpr
stress0(j,i,2,k) = s22v(k) +sarc2v*dtpr
stress0(j,i,4,k) = s33v(k) +sarc3v*dtpr ! Reset s33 for accretion M. Behn 04/19/05
endif

endif
!!!!!!!!!!!!!!!!!!!!!ends!!!!!!!!!!!!!!!!!!!
                end do
                irheol_fl (j,i) = -1
            endif
        endif

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
        endif

!  ACCUMULATED PLASTIC STRAIN
! Average the strain for pair of the triangles
! Note that area (n,it) is inverse of double area !!!!!
           if (sII_plas .lt. sII_visc) then
               aps(j,i) = aps(j,i) &
               + 0.5*( depl(1)*area(j,i,2)+depl(2)*area(j,i,1) ) / (area(j,i,1)+area(j,i,2)) &
               + 0.5*( depl(3)*area(j,i,4)+depl(4)*area(j,i,3) ) / (area(j,i,3)+area(j,i,4))
               if( aps(j,i) .lt. 0. ) aps(j,i) = 0.
           end if

! LINEAR HEALING OF THE PLASTIC STRAIN
           if (tau_heal .ne. 0.) &
               aps (j,i) = aps (j,i)/(1.+dt/tau_heal)
        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain(j,i,2) = strain(j,i,2) + 0.25*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain(j,i,3) = strain(j,i,3) + 0.25*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))

3 continue
!$OMP end do
!$OMP end parallel

if (inj_count.eq.inj_period+ext_period) inj_count = 0

devmax = max(devmax, curr_devmax)
dvmax = max(dvmax, curr_dvmax)


return
end
