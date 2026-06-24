module global_inputs

  !=========================
  ! Flags and global inputs
  ! Module added by A.JANIN 27.02.2026
  !=========================

  implicit none

  !=========================
  ! Global input arrays
  !=========================
  integer, allocatable :: i_kode(:)
  real*4, allocatable  :: i_ss(:)
  real*4, allocatable  :: i_ds(:)
  real*4, allocatable  :: i_ts(:)
  integer, allocatable :: i_fcode(:)
  real*4, allocatable  :: i_sdrop(:)
  real*4, allocatable  :: i_rhoLitho(:)
  real*4, allocatable  :: i_rhoFluid(:)
  real*4, allocatable  :: i_cohes(:)
  real*4, allocatable  :: i_disfric(:)
  integer, allocatable :: i_dyndike(:)
  real*4, allocatable  :: i_rhoMagma(:)
  real*4, allocatable  :: i_Hl(:)
  real*4, allocatable  :: i_DPm0(:)

  !=========================
  ! Generic solver parameters
  !=========================
  real(8), parameter :: EPS = 1000.0d0*epsilon(1.0d0)   ! Solver tolerance for the matrix inversion (in double precision)

  !=========================
  ! Friction solver parameters
  !=========================
  integer*4 :: i_maxiter
  real*4    :: i_tolsolver

  !=========================
  ! Global parameters
  !=========================
  real(8), parameter :: ggravi = 9.80665  ! average gravitational acceleration on Earth


contains

  subroutine allocate_global_inputs(ndis)
    implicit none
    integer, intent(in) :: ndis

    allocate(i_kode(ndis))
    allocate(i_ss(ndis))
    allocate(i_ds(ndis))
    allocate(i_ts(ndis))
    allocate(i_fcode(ndis))
    allocate(i_sdrop(ndis))
    allocate(i_rhoLitho(ndis))
    allocate(i_rhoFluid(ndis))
    allocate(i_cohes(ndis))
    allocate(i_disfric(ndis))
    allocate(i_dyndike(ndis))
    allocate(i_rhoMagma(ndis))
    allocate(i_Hl(ndis))
    allocate(i_DPm0(ndis))

  end subroutine allocate_global_inputs

  subroutine deallocate_global_inputs()
    implicit none

    deallocate(i_kode)
    deallocate(i_ss)
    deallocate(i_ds)
    deallocate(i_ts)
    deallocate(i_fcode)
    deallocate(i_sdrop)
    deallocate(i_rhoLitho)
    deallocate(i_rhoFluid)
    deallocate(i_cohes)
    deallocate(i_disfric)
    deallocate(i_dyndike)
    deallocate(i_rhoMagma)
    deallocate(i_Hl)
    deallocate(i_DPm0)

  end subroutine deallocate_global_inputs

end module global_inputs
