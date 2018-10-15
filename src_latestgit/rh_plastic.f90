
!------ Elasto-Plastic
!!!!!!!!!!!!!!!!!!!!!!!!!!!test
!subroutine mean(a,b,c,x,y,z)
!Implicit none
!real*8, intent(in) :: a,b,c
!real*8, intent(out) :: x,y,z
!x = a
!y = b
!z = c
!print *,'a,b,c =',a,b,c
!end subroutine
!!!!!!!!!!test ends!!!!!!!!!!!!
!!!!!!!!!!!

!!!!!!!!
subroutine plastic(bulkm,rmu,coh,phi,psi,depl,ipls,diss,hardn,s11,s22,s33,s12,s1a,s2a,s3a,de11,de22,de33,de12,&
     ten_off,ndim,irh_mark,x,y,coha,amca,anphia,phia,sphia,degrada,depla,fsa,e1a,e2a,ha,fta)
!subroutine plastic(bulkm,rmu,coh,phi,psi,depls,ipls,diss,hardn,s11,s22,s33,s12,de11,de22,de33,de12,&
!     ten_off,ndim,irh_mark)
implicit none

integer, intent(in) :: ndim, irh_mark,x,y
real*8, intent(in) :: bulkm, rmu, coh, phi, psi, hardn, de11, de22, de33, de12, ten_off
!real*8, intent(inout) :: s11, s22, s33, s12,s1,s2,s3
real*8, intent(inout) :: s11, s22, s33, s12!,s1a,s2a,s3a
real*8, intent(out) :: depl, diss,s1a,s2a,s3a,coha,anphia,amca,phia,sphia,degrada,depla,fsa,e1a,e2a,ha,fta
integer, intent(out) :: ipls
real*8, parameter :: pi = 3.14159265358979323846
real*8, parameter :: degrad = pi/180.
real*8, parameter :: c4d3 = 4./3.
real*8, parameter :: c2d3 = 2./3.
! press_add formaely was passed by a parameter. in my case it is always zero.
real*8, parameter :: press_add = 0. !4.e7
real*8 sphi, spsi, anphi, anpsi, amc, e1, e2, x1, ten_max, &
     s11i, s22i, s12i, s33i,s1,s2,s3, sdif, s0, rad, si, sii, psdif, &
     fs, alams, dep1, dep3, depm, cs2, si2, dc2, dss, sphii, anphii

!real*8 sphi, spsi, anphi, anpsi, amc, e1, e2, x1, ten_max, &
!     s11i, s22i, s12i, s33i, sdif, s0, rad, si, sii, psdif, &
!     fs, alams, dep1, dep3, depm, cs2, si2, dc2, dss, sphii, anphii
!integer icase
integer icase!,x,y

! ------------------------------
! Initialization section
! ------------------------------
depl = 0.
depla = 0.
diss = 0. 
ipls = 0 
      
sphi  = dsin(phi * degrad)
spsi  = dsin(psi * degrad)
!sphii = dcos(phi * degrad)
!anphii = sphii / (1. - sphi)
anphi = (1.+ sphi) / (1.- sphi)
anpsi = (1.+ spsi) / (1.- spsi)
amc   = 2.0 * coh * sqrt (anphi)
!amc =  2.0 * coh * anphii
e1    = bulkm + c4d3 * rmu
e2    = bulkm - c2d3 * rmu
x1    = (e1 - e2*anpsi + e1*anphi*anpsi - e2*anphi)
e1a = e1
e2a = e2
if (phi.eq. 0.) then
    ten_max=ten_off
else
    ten_max=min(ten_off,coh/(tan(phi*degrad)))
end if

! ---------------
! Running section
! ---------------

!---- get new trial stresses from old, assuming elastic increment
!---- add press (which is positive press = - (sxx+syy)*0.5, 
!---- which has 2 components: add pressure due to application of forces from the top 
!---- and subtract pressure of the fluid
s11i = s11 + (de22 + de33) *e2  + de11 *e1 - press_add
!open (unit = 1, file = "s11.txt")
!write (1,*) "Here are the s11 ", s11i
s22i = s22 + (de11 + de33) *e2  + de22 *e1 - press_add
!open (unit = 2, file = "s22.txt")
!write (2,*) "Here are the s22 ", s22i
s12i = s12 + de12 * 2.0 * rmu
s33i = s33 + (de11 + de22) *e2  + de33 *e1 - press_add
!open (unit = 3, file = "s33.txt")
!write (3,*) "Here are the s33 ", s33i
sdif = s11i - s22i
s0   = 0.5 * (s11i + s22i)
rad  = 0.5 * sqrt(sdif*sdif + 4.0 *s12i*s12i)
! principal stresses
si  = s0 - rad
sii = s0 + rad
psdif = si - sii
if (irh_mark.eq.1) then 
    s11 = s11i + press_add
    s22 = s22i + press_add
    s33 = s33i + press_add
    s12 = s12i
 return 
endif

!--------------------------------------------------------- 
!                         3D version 
!--------------------------------------------------------- 
if (ndim.eq.3) then
    !-- determine case ---
    if (s33i .gt. sii) then
        !- s33 is minor p.s. --
        icase = 3
        s1 = si
        s2 = sii
        s3 = s33i
    elseif (s33i .lt. si) then
        !- s33 is major p.s. --
        icase = 2
        s1 = s33i
        s2 = si
        s3 = sii
    else
        !- s33 is intermediate --
        icase = 1
        s1 = si
        s2 = s33i
        s3 = sii
    endif
endif

!------------------------------------------------------- 
!         2D version 
!------------------------------------------------------- 
if (ndim.eq.2) then 
    icase = 1
    s1 = si 
    s2 = s33i
    s3 = sii
    
endif 

!--------------------------------------------------------
! Check for tensional failure before the shear failure
!-------------------------------------------------------
 
!----- general tension failure


if (s1 .ge. ten_max) then
    ipls = -5
    goto 800
endif

!- uniaxial tension ... intermediate p.s. ---
if (s2 .ge. ten_max .and. ndim .eq.3) then
    ipls = -6
    s2 = ten_max
    s3 = ten_max
endif

!- partial failure (only if s3 is greater than ten_max) 
if (s3 .ge. ten_max) then
    s3 = ten_max 
    ipls = -7
endif

s1a = s1
s2a = s2
s3a = s3
coha = coh
amca = amc
anphia = anphi
phia = phi
sphia = sphi
degrada = degrad
fta = s3 - ten_off
!- check for shear yield (if fs<0 -> plastic flow)
fs = s1 - s3 * anphi + amc
fsa = fs
if (fs .lt. 0.0) then
    !-- yielding in shear ----
    if (icase .eq. 1) ipls = -2
    if (icase .eq. 2) ipls = -3
    if (icase .eq. 3) ipls = -4
    alams = fs/(x1+hardn)
    s1 = s1 - alams * (e1 - e2 * anpsi )
    s2 = s2 - alams * e2 * (1.0 - anpsi )
    s3 = s3 - alams * (e2 - e1 * anpsi )

    ! Increment of the plastic strain (2nd Invariant)
    dep1 = alams
    dep3 = -alams*anpsi

    ! FOR 2D caculations
    !dont touch this part
!!!!!!!!!origial depl!!!!!!!!!!!!!
    depm = 0.5*(dep1+dep3)
    depl = 0.5*abs(dep1-dep3)
    depla = depl
!!!!!!!!!!orignial ends!!!!!!!!!
    !depl = 0.5*abs(dep1+dep3)
    ! Dissipation rate
    diss = s1*dep1+s3*dep3
else
    !-- no failure at all (elastic behaviour)
    s11 = s11i + press_add
    s22 = s22i + press_add
    s33 = s33i + press_add
    s12 = s12i
    depla = depl
    return
endif

!- general tension failure?
if (s1 .ge. ten_max) then
    ipls = -5
    goto 800
endif

!- uniaxial tension ... intermediate p.s. ---
if (s2 .ge. ten_max .and.ndim.eq.3) then
    ipls = -6
    s2 = ten_max
    s3 = ten_max
    goto 205
endif

!- uniaxial tension ... minor p.s. ---
if (s3 .ge. ten_max) then
    ipls = -7
    s3 = ten_max 
endif

!- direction cosines
205 continue
if ( psdif .eq. 0. ) then
    cs2 = 1.
    si2 = 0.
else
    cs2 = sdif / psdif
    si2 = 2.0 * s12i / psdif
endif
s1a = s1
s2a = s2
s3a = s3
coha = coh
amca = amc
anphia = anphi
phia = phi
sphia = sphi
degrada = degrad
fta = s3 - ten_off
!print *, 's1a =', s1, 's2a =', s2,'s3a =', s3
!print *, 'si =', s1, s2, s3
!- resolve back to global axes
goto (210,220,230), icase

210 continue
dc2 = (s1-s3) * cs2
dss = s1 + s3
s11 = 0.5 * (dss + dc2)
s22 = 0.5 * (dss - dc2)
s12 = 0.5 * (s1 - s3) * si2
s33 = s2
goto 240

220 continue
dc2 = (s2-s3) * cs2
dss = s2 + s3
s11 = 0.5 * (dss + dc2)
s22 = 0.5 * (dss - dc2)
s12 = 0.5 * (s2 - s3) * si2
s33 = s1
goto 240

230 continue
dc2 = (s1-s2) * cs2
dss = s1 + s2
s11 = 0.5 * (dss + dc2)
s22 = 0.5 * (dss - dc2)
s12 = 0.5 * (s1 - s2) * si2
s33 = s3

240 continue

s11 = s11 + press_add
s22 = s22 + press_add
s33 = s33 + press_add

return

!-- set stresses to plastic apex ---
800   continue
s11        = ten_max
s22        = ten_max
s12        = 0.0
s33        = ten_max

s11 = s11 + press_add
s22 = s22 + press_add
s33 = s33 + press_add
return
end




!==================================================================
! Prepare plastic properties depending on softening, weighted by phase ratio

subroutine pre_plast (i,j,coh,phi,psi,hardn)
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

pls_curr = aps(j,i)

phi = 0
coh = 0
psi = 0

! Strain-Hardening
hardn = 0

do iph = 1, nphase
    if(phase_ratio(iph,j,i) .lt. 0.01) cycle

    if(pls_curr < plstrain1(iph)) then
        ! no weakening yet
        f = fric1(iph)
        c = cohesion1(iph)
        d = dilat1(iph)
        h = 0
    else if (pls_curr < plstrain2(iph)) then
        ! Find current properties from linear interpolation
        dpl = (pls_curr - plstrain1(iph)) / (plstrain2(iph) - plstrain1(iph))
        f =  fric1(iph) + (fric2(iph) - fric1(iph)) * dpl
        d = dilat1(iph) + (dilat2(iph) - dilat1(iph)) * dpl
        c = cohesion1(iph) + (cohesion2(iph) - cohesion1(iph)) * dpl
        h = (cohesion2(iph) - cohesion1(iph)) / (plstrain2(iph) - plstrain1(iph))
    else
        ! saturated weakening
        f = fric2(iph)
        c = cohesion2(iph)
        d = dilat2(iph)
        h = 0
    endif

    ! using harmonic mean on friction and cohesion
    ! using arithmatic mean on dilation and hardening
    phi = phi + phase_ratio(iph,j,i) / f
    coh = coh + phase_ratio(iph,j,i) / c
    psi = psi + phase_ratio(iph,j,i) * d
    hardn = hardn + phase_ratio(iph,j,i) * h
    ha = h
enddo

phi = 1 / phi
coh = 1 / coh
!if( mod(nloop, 10000).eq.0 .and.i.eq.28 .and.j .eq.5 ) then
!open (unit = 32, file = "cohm1.txt",position='append')
!write (32,*) coh, phi
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.27 .and.j .eq.5 ) then
!open (unit = 33, file = "cohl1.txt",position='append')
!write (33,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.29 .and.j .eq.5 ) then
!open (unit = 34, file = "cohr1.txt",position='append')
!write (34,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.28 .and.j .eq.4 ) then
!open (unit = 35, file = "coht1.txt",position='append')
!write (35,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.28 .and.j .eq.6 ) then
!open (unit = 36, file = "cohb1.txt",position='append')
!write (36,*) coh
!end if
!
!if( mod(nloop, 10000).eq.0 .and.i.eq.62 .and.j .eq.9 ) then
!open (unit = 37, file = "cohm2.txt",position='append')
!write (37,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.61 .and.j .eq.9 ) then
!open (unit = 38, file = "cohl2.txt",position='append')
!write (38,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.63 .and.j .eq.9 ) then
!open (unit = 39, file = "cohr2.txt",position='append')
!write (39,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.62 .and.j .eq.8 ) then
!open (unit = 40, file = "coht2.txt",position='append')
!write (40,*) coh
!end if
!if( mod(nloop, 10000).eq.0 .and.i.eq.62 .and.j .eq.10 ) then
!open (unit = 41, file = "cohb2.txt",position='append')
!write (41,*) coh
!end if
return
end subroutine pre_plast