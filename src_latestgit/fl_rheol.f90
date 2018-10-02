
!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  

subroutine fl_rheol
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension depl(4)
dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4),ps1(4),ps2(4),ps3(4),coha(4),amca(4),anphia(4)
dimension phia(4),sphia(4),degrada(4),depla(4),fsa(4),e1a(4),e2a(4),ha(4)
!dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
!real*8 bear1, bear2, bear3
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
if(irh.ge.11) call init_visc

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
    if(ny_inject.gt.0.and.j.le.nelem_inject) then
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



!        if (i .eq. 26 .and. j .eq. 5) then
!        Print *, i,j,coh,phi,psi




        ! Re-evaluate viscosity
        if (irh.eq.7 .or. irh.eq.12) then
!if (irh.eq.7 .or. irh.eq.12) then
            if( mod(nloop,ifreq_visc).eq.0 .OR. ireset.eq.1 ) visn(j,i) = Eff_visc(j,i)
        !            if (ny_inject.gt.0.and.i.eq.iinj) visn(j,i) = v_min
        endif
        !endif
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
!                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
!                    s11p(k),s22p(k),s33p(k),s12p(k),&
!                    de11,de22,de33,de12,ten_off,ndim,irh_mark)

!                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
!                    s11p(k),s22p(k),s33p(k),s12p(k),bear1, bear2, bear3,&
!                    de11,de22,de33,de12,ten_off,ndim,irh_mark,i,j)

                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                    s11p(k),s22p(k),s33p(k),s12p(k),ps1(k),ps2(k),ps3(k),&
                    de11,de22,de33,de12,ten_off,ndim,irh_mark,i,j,coha(k),&
                    amca(k),anphia(k),phia(k),sphia(k),degrada(k),depla(k),&
                    fsa(k),e1a(k),e2a(k),ha(k))

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
!                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
!                        s11p(k),s22p(k),s33p(k),s12p(k),&
!                        de11,de22,de33,de12,ten_off,ndim,irh_mark)

!                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
!                        s11p(k),s22p(k),s33p(k),s12p(k),bear1, bear2, bear3,&
!                        de11,de22,de33,de12,ten_off,ndim,irh_mark,i,j)
                        !print *, 'bear1 =', bear1, 'bear2 =', bear2,'bear3 =', bear3
!                        ps1(j,i,k) = bear1
!                        ps2(j,i,k) = bear2
!                        ps3(j,i,k) = bear3
                        !if( mod(nloop, 5000).eq.0 ) then
                        !print *, 'bear1r =', bear1, 'bear2r =', bear2,'bear3r =', bear3
                        !endif
                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                            s11p(k),s22p(k),s33p(k),s12p(k),ps1(k),ps2(k),ps3(k),&
                            de11,de22,de33,de12,ten_off,ndim,irh_mark,i,j,coha(k),&
                            amca(k),anphia(k),phia(k),sphia(k),degrada(k),depla(k),&
                            fsa(k),e1a(k),e2a(k),ha(k))
                    !!!!!!!!!!!!!!!assign priciple stresses for debugging!!!!!!!!!!!!!!!!
                        pss1(j,i,k) = ps1(k)
                        pss2(j,i,k) = ps2(k)
                        pss3(j,i,k) = ps3(k)
                        cohe(j,i,k) = coha(k)
                        amce(j,i,k) = amca(k)
                        anphie(j,i,k) = anphia(k)
                        phie(j,i,k) = phia(k)
                        sphie(j,i,k) = sphia(k)
                        degrade(j,i,k) = degrada(k)
                        deple(j,i,k) = depla(k)
                        fse(j,i,k) = fsa(k)
                        e1e(j,i,k) = e1a(k)
                        e2e(j,i,k) = e2a(k)
                        he(j,i,k) = ha(k)
!                        if( mod(nloop, 10000).eq.0 ) then
!                        print *, 'cohe =', cohe(5,28,1)
!                        endif
                        !print *, 'coh =',coh
!                        pss1(j,i,k) = ps1
!                        pss2(j,i,k) = ps2
!                        pss3(j,i,k) = ps3
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!debug ends!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                        de11,de22,de33,de12,dv,&
                        ndim,dt,curr_devmax,curr_dvmax)
                else ! use previously defined rheology
                    if( irheol_fl(j,i) .eq. 1 ) then
!                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
!                            s11p(k),s22p(k),s33p(k),s12p(k),&
!                            de11,de22,de33,de12,ten_off,ndim,irh_mark)

!                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
!                            s11p(k),s22p(k),s33p(k),s12p(k),bear1, bear2, bear3,&
!                            de11,de22,de33,de12,ten_off,ndim,irh_mark,i,j)
                         call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                             s11p(k),s22p(k),s33p(k),s12p(k),ps1(k),ps2(k),ps3(k),&
                             de11,de22,de33,de12,ten_off,ndim,irh_mark,i,j,coha(k),&
                             amca(k),anphia(k),phia(k),sphia(k),degrada(k),depla(k),&
                             fsa(k),e1a(k),e2a(k),ha(k))
!print *, 'bear1 =', bear1, 'bear3 =', bear2,'bear3 =', bear3
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
                    !print *,'ps1 =' ps1(k)
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

        !print *,'12'
        if (irh.eq.6 .or. irh.ge.11) then
            !  ACCUMULATED PLASTIC STRAIN
            ! Average the strain for pair of the triangles
            ! Note that area (n,it) is inverse of double area !!!!!
            dps = (0.5*( depl(1)*area(j,i,2)+depl(2)*area(j,i,1) ) / (area(j,i,1)+area(j,i,2))) &
                + (0.5*( depl(3)*area(j,i,4)+depl(4)*area(j,i,3) ) / (area(j,i,3)+area(j,i,4)))
            aps(j,i) = aps(j,i) + dps 
            if( aps(j,i) .lt. 0. ) aps(j,i) = 0.
            !print *, aps(5,28)
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

!!!!!!!!!!!!1st element of interests!!!!!!!!!!!!!
sec_year = 3.1558e+7
T_in_ma = time/sec_year/1.e6
ps1_avem1 = 0.25 * (pss1(5,28,1)+pss1(5,28,2)+pss1(5,28,3)+pss1(5,28,4))
!print*,'ps1ave =',ps1_avem1
ps2_avem1 = 0.25 * (pss2(5,28,1)+pss2(5,28,2)+pss2(5,28,3)+pss2(5,28,4))
ps3_avem1 = 0.25 * (pss3(5,28,1)+pss3(5,28,2)+pss3(5,28,3)+pss3(5,28,4))

!ps1_avel1 = 0.25 * (pss1(5,27,1)+pss1(5,27,2)+pss1(5,27,3)+pss1(5,27,4))
!ps2_avel1 = 0.25 * (pss2(5,27,1)+pss2(5,27,2)+pss2(5,27,3)+pss2(5,27,4))
!ps3_avel1 = 0.25 * (pss3(5,27,1)+pss3(5,27,2)+pss3(5,27,3)+pss3(5,27,4))

!ps1_aver1 = 0.25 * (pss1(5,29,1)+pss1(5,29,2)+pss1(5,29,3)+pss1(5,29,4))
!ps2_aver1 = 0.25 * (pss2(5,29,1)+pss2(5,29,2)+pss2(5,29,3)+pss2(5,29,4))
!ps3_aver1 = 0.25 * (pss3(5,29,1)+pss3(5,29,2)+pss3(5,29,3)+pss3(5,29,4))

ps1_avet1 = 0.25 * (pss1(5,27,1)+pss1(5,27,2)+pss1(5,27,3)+pss1(5,27,4))
ps2_avet1 = 0.25 * (pss2(5,27,1)+pss2(5,27,2)+pss2(5,27,3)+pss2(5,27,4))
ps3_avet1 = 0.25 * (pss3(5,27,1)+pss3(5,27,2)+pss3(5,27,3)+pss3(5,27,4))

!ps1_aveb1 = 0.25 * (pss1(6,28,1)+pss1(6,28,2)+pss1(6,28,3)+pss1(6,28,4))
!ps2_aveb1 = 0.25 * (pss2(6,28,1)+pss2(6,28,2)+pss2(6,28,3)+pss2(6,28,4))
!ps3_aveb1 = 0.25 * (pss3(6,28,1)+pss3(6,28,2)+pss3(6,28,3)+pss3(6,28,4))

coh_avem1 = 0.25 * (cohe(5,28,1)+cohe(5,28,2)+cohe(5,28,3)+cohe(5,28,4))
!coh_avel1 = 0.25 * (cohe(5,27,1)+cohe(5,27,2)+cohe(5,27,3)+cohe(5,27,4))
!coh_aver1 = 0.25 * (cohe(5,29,1)+cohe(5,29,2)+cohe(5,29,3)+cohe(5,29,4))
coh_avet1 = 0.25 * (cohe(5,27,1)+cohe(5,27,2)+cohe(5,27,3)+cohe(5,27,4))
!coh_aveb1 = 0.25 * (cohe(6,28,1)+cohe(6,28,2)+cohe(6,28,3)+cohe(6,28,4))

depl_avem1 = 0.25 * (deple(5,28,1)+deple(5,28,2)+deple(5,28,3)+deple(5,28,4))
depl_avet1 = 0.25 * (deple(5,27,1)+deple(5,27,2)+deple(5,27,3)+deple(5,27,4))
acum_plsm1 =  aps(5,28)
acum_plst1 =  aps(5,27)
fse_avem1 = 0.25 * (fse(5,28,1)+fse(5,28,2)+fse(5,28,3)+fse(5,28,4))
fse_avet1 = 0.25 * (fse(5,27,1)+fse(5,27,2)+fse(5,27,3)+fse(5,27,4))
e1e_avem1 = 0.25 * (e1e(5,28,1)+e1e(5,28,2)+e1e(5,28,3)+e1e(5,28,4))
e1e_avet1 = 0.25 * (e1e(5,27,1)+e1e(5,27,2)+e1e(5,27,3)+e1e(5,27,4))
e2e_avem1 = 0.25 * (e2e(5,28,1)+e2e(5,28,2)+e2e(5,28,3)+e2e(5,28,4))
e2e_avet1 = 0.25 * (e2e(5,27,1)+e2e(5,27,2)+e2e(5,27,3)+e2e(5,27,4))
he_avem1 = 0.25 * (he(5,28,1)+he(5,28,2)+he(5,28,3)+he(5,28,4))
he_avet1 = 0.25 * (he(5,27,1)+he(5,27,2)+he(5,27,3)+he(5,27,4))
amc_avem1 = 0.25 * (amce(5,28,1)+amce(5,28,2)+amce(5,28,3)+amce(5,28,4))
amc_avel1 = 0.25 * (amce(5,27,1)+amce(5,27,2)+amce(5,27,3)+amce(5,27,4))
!amc_aver1 = 0.25 * (amce(5,29,1)+amce(5,29,2)+amce(5,29,3)+amce(5,29,4))
!amc_avet1 = 0.25 * (amce(4,28,1)+amce(4,28,2)+amce(4,28,3)+amce(4,28,4))
!amc_aveb1 = 0.25 * (amce(6,28,1)+amce(6,28,2)+amce(6,28,3)+amce(6,28,4))
!
!
anphi_avem1 = 0.25 * (anphie(5,28,1)+anphie(5,28,2)+anphie(5,28,3)+anphie(5,28,4))
anphi_avel1 = 0.25 * (anphie(5,27,1)+anphie(5,27,2)+anphie(5,27,3)+anphie(5,27,4))
!anphi_aver1 = 0.25 * (anphie(5,29,1)+anphie(5,29,2)+anphie(5,29,3)+anphie(5,29,4))
!anphi_avet1 = 0.25 * (anphie(4,28,1)+anphie(4,28,2)+anphie(4,28,3)+anphie(4,28,4))
!anphi_aveb1 = 0.25 * (anphie(6,28,1)+anphie(6,28,2)+anphie(6,28,3)+anphie(6,28,4))
!
!phi_avem1 = 0.25 * (phie(5,28,1)+phie(5,28,2)+phie(5,28,3)+phie(5,28,4))
!phi_avel1 = 0.25 * (phie(5,27,1)+phie(5,27,2)+phie(5,27,3)+phie(5,27,4))
!phi_aver1 = 0.25 * (phie(5,29,1)+phie(5,29,2)+phie(5,29,3)+phie(5,29,4))
!phi_avet1 = 0.25 * (phie(4,28,1)+phie(4,28,2)+phie(4,28,3)+phie(4,28,4))
!phi_aveb1 = 0.25 * (phie(6,28,1)+phie(6,28,2)+phie(6,28,3)+phie(6,28,4))
!
!sphi_avem1 = 0.25 * (sphie(5,28,1)+sphie(5,28,2)+sphie(5,28,3)+sphie(5,28,4))
!sphi_avet1 = 0.25 * (sphie(4,28,1)+sphie(4,28,2)+sphie(4,28,3)+sphie(4,28,4))
!degrad_avem1 = 0.25 * (degrade(5,28,1)+degrade(5,28,2)+degrade(5,28,3)+degrade(5,28,4))
!degrad_avet1 = 0.25 * (degrade(4,28,1)+degrade(4,28,2)+degrade(4,28,3)+degrade(4,28,4))
!
!!!!!!!!!!!!!!!!!!2nd element of interests!!!!!!!!!!!!!!!
!ps1_avem2 = 0.25 * (pss1(9,62,1)+pss1(9,62,2)+pss1(9,62,3)+pss1(9,62,4))
!ps2_avem2 = 0.25 * (pss2(9,62,1)+pss2(6,62,2)+pss2(9,62,3)+pss2(9,62,4))
!ps3_avem2 = 0.25 * (pss3(9,62,1)+pss3(6,62,2)+pss3(9,62,3)+pss3(9,62,4))
!
!ps1_avel2 = 0.25 * (pss1(9,61,1)+pss1(9,61,2)+pss1(9,61,3)+pss1(9,61,4))
!ps2_avel2 = 0.25 * (pss2(9,61,1)+pss2(9,61,2)+pss2(9,61,3)+pss2(9,61,4))
!ps3_avel2 = 0.25 * (pss3(9,61,1)+pss3(9,61,2)+pss3(9,61,3)+pss3(9,61,4))
!
!ps1_aver2 = 0.25 * (pss1(9,63,1)+pss1(9,63,2)+pss1(9,63,3)+pss1(9,63,4))
!ps2_aver2 = 0.25 * (pss2(9,63,1)+pss2(9,63,2)+pss2(9,63,3)+pss2(9,63,4))
!ps3_aver2 = 0.25 * (pss3(9,63,1)+pss3(9,63,2)+pss3(9,63,3)+pss3(9,63,4))
!
!ps1_avet2 = 0.25 * (pss1(8,62,1)+pss1(8,62,2)+pss1(8,62,3)+pss1(8,62,4))
!ps2_avet2 = 0.25 * (pss2(8,62,1)+pss2(8,62,2)+pss2(8,62,3)+pss2(8,62,4))
!ps3_avet2 = 0.25 * (pss3(8,62,1)+pss3(8,62,2)+pss3(8,62,3)+pss3(8,62,4))
!
!ps1_aveb2 = 0.25 * (pss1(10,62,1)+pss1(10,62,2)+pss1(10,62,3)+pss1(10,62,4))
!ps2_aveb2 = 0.25 * (pss2(10,62,1)+pss2(10,62,2)+pss2(10,62,3)+pss2(10,62,4))
!ps3_aveb2 = 0.25 * (pss3(10,62,1)+pss3(10,62,2)+pss3(10,62,3)+pss3(10,62,4))
!!
!if( (T_in_ma.ge.0.31) .and. (i .eq. 28) .and. (j .eq. 5) ) then
!open (unit = 1, file = "ps1m1.txt",position='append')
!open (unit = 2, file = "ps2m1.txt",position='append')
!open (unit = 3, file = "ps3m1.txt",position='append')
!write (1,*) ps1_avem1
!write (2,*) ps2_avem1
!write (3,*) ps3_avem1
!
!open (unit = 32, file = "cohm1.txt",position='append')
!write (32,*) coh_avem1
!
!open (unit = 56, file = "deplsm1.txt",position='append')
!write (56,*) depl_avem1
!
!open (unit = 58, file = "acuplsm1.txt",position='append')
!write (58,*) acum_plsm1
!
!open (unit = 60, file = "fsm1.txt",position='append')
!write (60,*) fse_avem1
!
!open (unit = 62, file = "e1m1.txt",position='append')
!write (62,*) e1e_avem1
!
!
!open (unit = 64, file = "e2m1.txt",position='append')
!write (64,*) e2e_avem1
!
!open (unit = 66, file = "hm1.txt",position='append')
!write (66,*) he_avem1
!
!open (unit = 37, file = "amcm1.txt",position='append')
!write (37,*) amc_avem1
!
!open (unit = 42, file = "anphim1.txt",position='append')
!write (42,*) anphi_avem1
!
!open (unit = 16, file = "timeinmy.txt",position='append')
!write (16,*) (time/sec_year/1.e6)
!end if
!
!if( (T_in_ma.ge.0.31) .and. (i .eq. 27) .and. (j .eq. 5) ) then
!!open (unit = 4, file = "ps1l1.txt",position='append')
!!open (unit = 5, file = "ps2l1.txt",position='append')
!!open (unit = 6, file = "ps3l1.txt",position='append')
!!write (4,*) ps1_avel1
!!write (5,*) ps2_avel1
!!write (6,*) ps3_avel1
!!
!!open (unit = 7, file = "ps1r1.txt",position='append')
!!open (unit = 8, file = "ps2r1.txt",position='append')
!!open (unit = 9, file = "ps3r1.txt",position='append')
!!write (7,*) ps1_aver1
!!write (8,*) ps2_aver1
!!write (9,*) ps3_aver1
!
!open (unit = 10, file = "ps1t1.txt",position='append')
!open (unit = 11, file = "ps2t1.txt",position='append')
!open (unit = 12, file = "ps3t1.txt",position='append')
!write (10,*) ps1_avet1
!write (11,*) ps2_avet1
!write (12,*) ps3_avet1
!
!!open (unit = 13, file = "ps1b1.txt",position='append')
!!open (unit = 14, file = "ps2b1.txt",position='append')
!!open (unit = 15, file = "ps3b1.txt",position='append')
!!write (13,*) ps1_aveb1
!!write (14,*) ps2_aveb1
!!write (15,*) ps3_aveb1
!
!
!!open (unit = 33, file = "cohl1.txt",position='append')
!!write (33,*) coh_avel1
!!open (unit = 34, file = "cohr1.txt",position='append')
!!write (34,*) coh_aver1
!open (unit = 35, file = "coht1.txt",position='append')
!write (35,*) coh_avet1
!!open (unit = 36, file = "cohb1.txt",position='append')
!!write (36,*) coh_aveb1
!
!
!!open (unit = 38, file = "amcl1.txt",position='append')
!!write (38,*) amc_avel1
!!open (unit = 39, file = "amcr1.txt",position='append')
!!write (39,*) amc_aver1
!open (unit = 40, file = "amct1.txt",position='append')
!write (40,*) amc_avet1
!!open (unit = 41, file = "amcb1.txt",position='append')
!!write (41,*) amc_aveb1
!!
!
!!open (unit = 43, file = "anphil1.txt",position='append')
!!write (43,*) anphi_avel1
!!open (unit = 44, file = "anphir1.txt",position='append')
!!write (44,*) anphi_aver1
!open (unit = 45, file = "anphit1.txt",position='append')
!write (45,*) anphi_avet1
!!open (unit = 46, file = "anphib1.txt",position='append')
!!write (46,*) anphi_aveb1
!!
!!open (unit = 47, file = "phim1.txt",position='append')
!!write (47,*) phi_avem1
!!open (unit = 48, file = "phil1.txt",position='append')
!!write (48,*) phi_avel1
!!open (unit = 49, file = "phir1.txt",position='append')
!!write (49,*) phi_aver1
!!open (unit = 50, file = "phit1.txt",position='append')
!!write (50,*) phi_avet1
!!open (unit = 51, file = "phib1.txt",position='append')
!!write (51,*) phi_aveb1
!!
!!
!!open (unit = 52, file = "sphim1.txt",position='append')
!!write (52,*) sphi_avem1
!!
!!open (unit = 53, file = "sphit1.txt",position='append')
!!write (53,*) sphi_avet1
!!
!!open (unit = 54, file = "degradm1.txt",position='append')
!!write (54,*) degrad_avem1
!!
!!open (unit = 55, file = "degradt1.txt",position='append')
!!write (55,*) degrad_avet1
!
!open (unit = 57, file = "deplst1.txt",position='append')
!write (57,*) depl_avet1
!
!
!open (unit = 59, file = "acuplst1.txt",position='append')
!write (59,*) acum_plst1
!
!
!
!
!
!open (unit = 61, file = "fst1.txt",position='append')
!write (61,*) fse_avet1
!
!
!open (unit = 63, file = "e1t1.txt",position='append')
!write (63,*) e1e_avet1
!
!
!open (unit = 65, file = "e2t1.txt",position='append')
!write (65,*) e2e_avet1
!
!
!open (unit = 67, file = "ht1.txt",position='append')
!write (67,*) he_avet1
!
!end if
!
!!!
!!!
!!!
!!open (unit = 17, file = "ps1m2.txt",position='append')
!!open (unit = 18, file = "ps2m2.txt",position='append')
!!open (unit = 19, file = "ps3m2.txt",position='append')
!!write (17,*) ps1_avem2
!!write (18,*) ps2_avem2
!!write (19,*) ps3_avem2
!!
!!open (unit = 20, file = "ps1l2.txt",position='append')
!!open (unit = 21, file = "ps2l2.txt",position='append')
!!open (unit = 22, file = "ps3l2.txt",position='append')
!!write (20,*) ps1_avel2
!!write (21,*) ps2_avel2
!!write (22,*) ps3_avel2
!!
!!open (unit = 23, file = "ps1r2.txt",position='append')
!!open (unit = 24, file = "ps2r2.txt",position='append')
!!open (unit = 25, file = "ps3r2.txt",position='append')
!!write (23,*) ps1_aver2
!!write (24,*) ps2_aver2
!!write (25,*) ps3_aver2
!!
!!open (unit = 26, file = "ps1t2.txt",position='append')
!!open (unit = 27, file = "ps2t2.txt",position='append')
!!open (unit = 28, file = "ps3t2.txt",position='append')
!!write (26,*) ps1_avet2
!!write (27,*) ps2_avet2
!!write (28,*) ps3_avet2
!!
!!open (unit = 29, file = "ps1b2.txt",position='append')
!!open (unit = 30, file = "ps2b2.txt",position='append')
!!open (unit = 31, file = "ps3b2.txt",position='append')
!!write (29,*) ps1_aveb2
!!write (30,*) ps2_aveb2
!!write (31,*) ps3_aveb2
!!!

!!!!!!!!!!!debug ends!!!!!!!!!!!!!!


        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain(j,i,2) = strain(j,i,2) + 0.25*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain(j,i,3) = strain(j,i,3) + 0.25*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))

3 continue
!$OMP end do
!$OMP end parallel


devmax = max(devmax, curr_devmax)
dvmax = max(dvmax, curr_dvmax)

return
end
