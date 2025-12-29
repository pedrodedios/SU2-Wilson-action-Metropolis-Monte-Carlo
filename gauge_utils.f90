!====================================================================================================
!
! Module: gauge_utils
!
! Description:
!   This module provides utility routines for SU(2) lattice gauge theory simulations.
!   It contains basic lattice indexing tools, random SU(2) matrix generators, and
!   low-level linear algebra operations for 2×2 complex matrices.
!
!   The utilities implemented here are used throughout the Monte Carlo simulation
!   and are designed to be lightweight, reusable, and numerically stable.
!
!   The module includes:
!     - Index shift functions with periodic boundary conditions
!     - Mapping of 4D lattice coordinates to a 1D site index
!     - Generation of random SU(2) matrices using the Haar measure
!     - Generation of small SU(2) rotations near the identity
!     - SU(2) matrix algebra (multiplication, dagger, determinant, trace)
!
!   All SU(2) matrices are represented as 2×2 complex matrices in double precision.
!
! Usage:
!   This module must be compiled before executing the main program, wilson_gauge_main.f90.
!
!====================================================================================================






module gauge_utils




implicit none
integer, parameter :: dp = kind(1.0d0)
  
contains





!shifts of indices in up and down directions (periodic boundary conditions)
integer function shiftup(num, NS)
    implicit none
    integer, intent(in) :: num, NS
    shiftup = num + 1
    if (num == NS) shiftup = 1
end function shiftup

integer function shiftdown(num, NS)
    implicit none
    integer, intent(in) :: num, NS
    shiftdown = num - 1
    if (num == 1) shiftdown = NS
end function shiftdown






!Master 1D index 
integer function toIndex(x, y, z, t, XS, YS, ZS, TS)
    implicit none
    integer, intent(in) :: x,y,z,t,XS,YS,ZS,TS
    toIndex = (x-1) + XS*((y-1) + YS*((z-1) + ZS*(t-1))) + 1
end function toIndex







subroutine SU2_HAAR(U)  ! hot start / U configuration proposal
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), intent(out) :: U(2,2)
    real(dp) :: a0, a1, a2, a3, norm

     call random_number(a0)
     call random_number(a1)
     call random_number(a2)
     call random_number(a3)

     ! Shift to (-0.5,0.5)
     a0 = a0 - 0.5_dp
     a1 = a1 - 0.5_dp
     a2 = a2 - 0.5_dp
     a3 = a3 - 0.5_dp

     norm = sqrt(a0*a0 + a1*a1 + a2*a2 + a3*a3)

     a0 = a0 / norm
     a1 = a1 / norm
     a2 = a2 / norm
     a3 = a3 / norm

     U(1,1) = cmplx(a0,  a3, dp)
     U(1,2) = cmplx(a2,  a1, dp)
     U(2,1) = cmplx(-a2, a1, dp)
     U(2,2) = cmplx(a0, -a3, dp)

end subroutine SU2_HAAR









! optianl method to construct U'= R U, where R is a small rotation
! near the identity controlled by the parameter eps.
subroutine SU2_RANDOM_NEAR_IDENTITY(R, eps)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  complex(dp), intent(out) :: R(2,2)
  real(dp), intent(in) :: eps

  real(dp) :: r1, r2, r3, norm, alpha
  real(dp) :: c, s

  ! random direction on S^2
  call random_number(r1)
  call random_number(r2)
  call random_number(r3)

  r1 = 2.0_dp*r1 - 1.0_dp
  r2 = 2.0_dp*r2 - 1.0_dp
  r3 = 2.0_dp*r3 - 1.0_dp

  norm = sqrt(r1*r1 + r2*r2 + r3*r3)
  r1 = r1 / norm
  r2 = r2 / norm
  r3 = r3 / norm

  ! small angle
  call random_number(alpha)
  alpha = eps * (2.0_dp*alpha - 1.0_dp)

  c = cos(alpha)
  s = sin(alpha)

  ! SU(2) matrix
  R(1,1) = cmplx(c,  s*r3, dp)
  R(1,2) = cmplx(s*r2, s*r1, dp)
  R(2,1) = cmplx(-s*r2, s*r1, dp)
  R(2,2) = cmplx(c, -s*r3, dp)

end subroutine SU2_RANDOM_NEAR_IDENTITY




! real trace of the SU(2) matrices
double precision function realTraceM(mat)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    complex(kind=dp), intent(in) :: mat(2,2)
    integer :: i

    realTraceM = 0.0d0
    do i = 1, 2
        realTraceM = realTraceM + real(mat(i,i), kind=dp)
    end do
end function realTraceM




! 2 matrices multiplication
subroutine matrix_mult(A, B, C)

    use iso_fortran_env, only: dp => real64
    implicit none
    complex(dp), intent(in)  :: A(2,2), B(2,2)
    complex(dp), intent(out) :: C(2,2)

    C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1)
    C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2)
    C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1)
    C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2)

end subroutine



! 3 matrices multiplication
pure subroutine matrix_mult3(A, B, C, R)
        
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), intent(in)  :: A(2,2), B(2,2), C(2,2)
    complex(dp), intent(out) :: R(2,2)

    complex(dp) :: t11, t12, t21, t22

    ! First T = B*C (kept in registers)
    t11 = B(1,1)*C(1,1) + B(1,2)*C(2,1)
    t12 = B(1,1)*C(1,2) + B(1,2)*C(2,2)
    t21 = B(2,1)*C(1,1) + B(2,2)*C(2,1)
    t22 = B(2,1)*C(1,2) + B(2,2)*C(2,2)

     ! Then R = A*T
     R(1,1) = A(1,1)*t11 + A(1,2)*t21
     R(1,2) = A(1,1)*t12 + A(1,2)*t22
     R(2,1) = A(2,1)*t11 + A(2,2)*t21
     R(2,2) = A(2,1)*t12 + A(2,2)*t22

end subroutine matrix_mult3




! 4 matrices multiplication
pure subroutine matrix_mult4(A, B, C, D, R)

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), intent(in)  :: A(2,2), B(2,2), C(2,2), D(2,2)
    complex(dp), intent(out) :: R(2,2)

    complex(dp) :: t11, t12, t21, t22
    complex(dp) :: s11, s12, s21, s22

    ! T = C*D
    t11 = C(1,1)*D(1,1) + C(1,2)*D(2,1)
    t12 = C(1,1)*D(1,2) + C(1,2)*D(2,2)
    t21 = C(2,1)*D(1,1) + C(2,2)*D(2,1)
    t22 = C(2,1)*D(1,2) + C(2,2)*D(2,2)

    ! S = B*T
    s11 = B(1,1)*t11 + B(1,2)*t21
    s12 = B(1,1)*t12 + B(1,2)*t22
    s21 = B(2,1)*t11 + B(2,2)*t21
    s22 = B(2,1)*t12 + B(2,2)*t22

    ! R = A*S
    R(1,1) = A(1,1)*s11 + A(1,2)*s21
    R(1,2) = A(1,1)*s12 + A(1,2)*s22
    R(2,1) = A(2,1)*s11 + A(2,2)*s21
    R(2,2) = A(2,1)*s12 + A(2,2)*s22

end subroutine matrix_mult4



!conjugate transpose matrix
pure subroutine dagger(A, B)

    implicit none
    complex(dp), intent(in) :: A(2,2)
    complex(dp), intent(out):: B(2,2)
    
    B(1,1) = conjg(A(1,1))
    B(1,2) = conjg(A(2,1))
    B(2,1) = conjg(A(1,2))
    B(2,2) = conjg(A(2,2))

end subroutine



function det(A) result(detA)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    complex(kind=dp), intent(in) :: A(2,2)
    complex(kind=dp) :: detA

    detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)
end function det









end module gauge_utils
