! -*- F90 -*-

module arrays
  implicit none

  ! fortran array pointer
  real*8, pointer, save :: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), balance(:,:,:), amass(:,:), rmass(:,:), &
area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:),pss1(:,:,:),pss2(:,:,:),pss3(:,:,:),cohe(:,:,:), &
amce(:,:,:), anphie(:,:,:),phie(:,:,:),degrade(:,:,:),sphie(:,:,:),deple(:,:,:), fse(:,:,:),&
e1e(:,:,:),e2e(:,:,:),he(:,:,:)

  integer, pointer, save :: ncod(:,:,:)

  ! temporary array
  real*8, pointer, save :: junk2(:,:)


contains

  subroutine allocate_arrays(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx

    allocate(cord(nz, nx, 2))
    allocate(temp(nz, nx))
    allocate(vel(nz, nx, 2))
    allocate(stress0(nz-1, nx-1, 4, 4))
    allocate(pss1(nz-1,nx-1,4))
    allocate(pss2(nz-1,nx-1,4))
    allocate(pss3(nz-1,nx-1,4))
    allocate(fse(nz-1,nx-1,4))
    allocate(e1e(nz-1,nx-1,4))
    allocate(e2e(nz-1,nx-1,4))
    allocate(he(nz-1,nx-1,4))
    allocate(cohe(nz-1,nx-1,4))
    allocate(amce(nz-1,nx-1,4))
    allocate(anphie(nz-1,nx-1,4))
    allocate(phie(nz-1,nx-1,4))
    allocate(anphie(nz-1,nx-1,4))
    allocate(degrade(nz-1,nx-1,4))
    allocate(sphie(nz-1,nx-1,4))
    allocate(deple(nz-1,nx-1,4))
    allocate(force(nz, nx, 2))
    allocate(balance(nz, nx, 2))
    allocate(amass(nz, nx))
    allocate(rmass(nz, nx))
    allocate(area(nz-1, nx-1, 4))
    allocate(dvol(nz-1, nx-1, 4))
    allocate(strain(nz-1, nx-1, 3))
    allocate(bc(nz, nx, 2))

    allocate(ncod(nz, nx, 2))

    allocate(junk2(nz, nx))

  end subroutine allocate_arrays

end module arrays
