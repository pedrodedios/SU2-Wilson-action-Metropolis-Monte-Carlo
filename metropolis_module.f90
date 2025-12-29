!====================================================================================================
!
! Module: metropolis_module
!
! Description:
!   This module implements the Metropolis Monte Carlo algorithm for SU(2) pure gauge theory
!   on a four-dimensional Euclidean lattice. The update procedure follows a local Metropolis
!   scheme with periodic boundary conditions and uses SU(2) link variables represented as
!   2×2 complex matrices.
!
!   The module performs:
!     - Hot-start initialization using the SU(2) Haar measure
!     - Local Metropolis updates of link variables
!     - Construction of staples for the Wilson gauge action
!     - Measurement of gauge-invariant observables (plaquette, action, susceptibility)
!     - Optional OpenMP parallelization for plaquette measurements
!
!
!   All SU(2) matrix operations (multiplication, dagger, trace, random generation)
!   are provided by the gauge_utils module.
!
! Usage:
!   This module is intended to be used by the main program wilson_gauge.f90.
!   Compile this module after gauge_utils.f90 and before the main program.
!
!   ***More involved observables can be further investigated using this module.***
!====================================================================================================







module metropolis




use gauge_utils    !math and SU(2) utilities
implicit none

contains




subroutine metropolis_gauge(beta, XS, YS, ZS, TS, eps, num_sweeps, thermalization, measuring_gap,u_action, u_plaq, u_susc)

    use gauge_utils
    implicit none



    ! parameters
    integer, intent(in) :: XS, YS, ZS,TS, num_sweeps, thermalization, measuring_gap
    real(8), intent(in) :: beta, eps
    integer :: VOL

    ! variables
    integer, allocatable :: uptable(:,:), downtable(:,:)
    integer :: site, x, y, z, t, mu, nu, sweep
    integer :: i, j, count
    integer :: site_minu, site_mu_min_nu, vertex
    double precision :: w, DELTA, S
    double precision :: deviation
    double precision, parameter :: PI = 4.0d0 * atan(1.0d0)
    double precision :: list_plaquette(1:50000), list_plaquette_sqr(1:50000)
    double precision :: plaquette_sum, plaquette_avg, plaquette_avg_square
    double precision :: plaquette_sum_square, trace_plaquette
    double precision :: temp_plaq, temp_sqr
    double precision :: mean_plaquette, mean_plaq_sqr
    double precision :: susceptibility, std_dev


    complex(kind=dp), allocatable :: configuration(:,:,:,:), prima(:,:,:,:)
    complex(kind=dp) :: identity(2,2), plaquette(2,2)
    complex(kind=dp) :: matrixproduct(2,2)
    complex(kind=dp) :: stap(2,2), temp_prima(2,2)
    complex(kind=dp) :: U1(2,2), U2(2,2), U3(2,2), U4(2,2)
    complex(kind=dp) :: U5(2,2), U6(2,2)
    complex(kind=dp) :: staple(2,2), matrix_plaquette(2,2)
    complex(kind=dp) :: local_contribution(2,2)
    complex(kind=dp) :: temp1(2,2), temp2(2,2)
    complex(kind=dp) :: up_staple(2,2), low_staple(2,2)

    integer :: u_action, u_plaq, u_susc


    VOL = XS * YS * ZS * TS

    allocate(uptable(VOL,4))
    allocate(downtable(VOL,4))
    allocate(configuration(VOL,4,2,2))
    allocate(prima(VOL,4,2,2))


    ! identity matrix
    do i = 1, 2
        do j = 1, 2
            if (i == j) then
                identity(i,j) = cmplx(1.0d0, 0.0d0)
            else
                identity(i,j) = cmplx(0.0d0, 0.0d0)
            end if
        end do
    end do


    ! initialize up and down tables of neighbour indices
    do x = 1, XS
        do y = 1, YS
            do z = 1, ZS
                do t = 1, TS
                    site = toIndex(x,y,z,t,XS,YS,ZS,TS)

                    uptable(site,1) = toIndex(shiftup(x,XS), y, z, t, XS, YS, ZS, TS)
                    uptable(site,2) = toIndex(x, shiftup(y,YS), z, t, XS, YS, ZS, TS)
                    uptable(site,3) = toIndex(x, y, shiftup(z,ZS), t, XS, YS, ZS, TS)
                    uptable(site,4) = toIndex(x, y, z, shiftup(t,TS), XS, YS, ZS, TS)

                    downtable(site,1) = toIndex(shiftdown(x,XS), y, z, t, XS, YS, ZS, TS)
                    downtable(site,2) = toIndex(x, shiftdown(y,YS), z, t, XS, YS, ZS, TS)
                    downtable(site,3) = toIndex(x, y, shiftdown(z,ZS), t, XS, YS, ZS, TS)
                    downtable(site,4) = toIndex(x, y, z, shiftdown(t,TS), XS, YS, ZS, TS)
                end do
            end do
        end do
    end do


    ! hot start
    do site = 1, VOL
        do mu = 1, 4
            call SU2_HAAR(configuration(site, mu, :, :))
        end do
    end do


    count = 0

    ! perform sweeps 
    do sweep = 1, num_sweeps

        do site = 1, VOL
            do mu = 1, 4

                call SU2_HAAR(temp_prima)             ! Random SU(2) proposal U'
                prima(site, mu, :, :) = temp_prima

                !call SU2_RANDOM_NEAR_IDENTITY(small_rotation, eps)    ** optional proposal algorithm, small rotation near identity **
                !call matrix_mult(small_rotation, configuration(site, mu, :, :), prima(site, mu, :, :)) ! propose new link U' = R*U 

                staple = (0.0d0, 0.0d0)

                do nu = 1, 4
                    if (mu /= nu) then

                        ! upper staple
                        U1 = configuration(uptable(site, mu), nu, :, :)                 ! Uν(x+μ)
                        call dagger(configuration(uptable(site, nu), mu, :, :), U2)     ! Uμ†(x+ν)
                        call dagger(configuration(site, nu, :, :), U3)                  ! Uν†(x)
                        call matrix_mult3(U1, U2, U3, up_staple)

                        ! lower staple
                        site_minu = downtable(site, nu)                                  ! x-ν
                        site_mu_min_nu = uptable(site_minu, mu)                          ! x+μ-ν     

                        call dagger(configuration(site_mu_min_nu, nu, :, :), U4)         ! Uν†(x+μ-ν)
                        call dagger(configuration(site_minu, mu, :, :), U5)              ! Uμ†(x-ν)
                        U6 = configuration(site_minu, nu, :, :)                          ! Uν(x-ν)
                        call matrix_mult3(U4, U5, U6, low_staple)

                        staple = staple + up_staple + low_staple

                    end if
                end do

                call matrix_mult(prima(site,mu,:,:), staple, temp1)
                call matrix_mult(configuration(site,mu,:,:), staple, temp2)

                DELTA = (-beta/2.0d0) * (realTraceM(temp1) - realTraceM(temp2))       ! compute change on the action 

                call RANDOM_NUMBER(w)

                if (DELTA <= 0.0d0 .or. log(w) < -DELTA) then                         ! Metropolis update
                    configuration(site, mu, :, :) = prima(site, mu, :, :)
                end if

            end do
        end do


        if (sweep > thermalization .and. mod(sweep, measuring_gap) == 0) then         ! Perform measuraments after thermalization



            count = count + 1                      ! # of measuraments
            S = 0.0d0                              ! action sum
            plaquette_sum = 0.0d0                  ! plaquette sum
            plaquette_sum_square = 0.0d0           ! plaquette**2 sum
 
            !Parallelized loops
            !$omp parallel do collapse(2) reduction(+:plaquette_sum,S) &                      
            !$omp private(U1,U2,U3,U4,matrix_plaquette,trace_plaquette,local_contribution)


            do vertex = 1, VOL
                do mu = 1, 4
                    do nu = mu+1, 4
!
                        U1 = configuration(vertex, mu, :, :)                          ! Uμ(x)
                        U2 = configuration(uptable(vertex, mu), nu, :, :)             ! Uν(x+μ)
                        call dagger(configuration(uptable(vertex,nu), mu, :, :), U3)  ! Uμ†(x+ν)
                        call dagger(configuration(vertex, nu, :, :), U4)              ! Uν†(x)

                        call matrix_mult4(U1, U2, U3, U4, matrix_plaquette)

                        local_contribution = identity - matrix_plaquette
                        S = S + (beta / 2.0d0) * realTraceM(local_contribution)       !action

                        trace_plaquette = 0.5d0 * realTraceM(matrix_plaquette)
                        plaquette_sum = plaquette_sum + trace_plaquette               ! plaquette 

                    end do
                end do
            end do

            !$omp end parallel do

            plaquette_avg = plaquette_sum / (6.0d0 * dble(VOL))                       !plaquette average
            plaquette_avg_square = plaquette_avg**2

            write(u_action,*) count, S                                                !file with Action data
            write(u_plaq,*) count, plaquette_avg, plaquette_avg_square                !file with plaquette data

            list_plaquette(count) = plaquette_avg
            list_plaquette_sqr(count) = plaquette_avg_square

        end if

    end do

    temp_plaq = 0.0d0
    temp_sqr  = 0.0d0

    do i = 1, count
        temp_plaq = temp_plaq + list_plaquette(i)
        temp_sqr  = temp_sqr  + list_plaquette_sqr(i)
    end do

    mean_plaquette = temp_plaq / count
    mean_plaq_sqr = temp_sqr / count

    susceptibility = (6.0d0 * dble(VOL)) * (mean_plaq_sqr - mean_plaquette**2)

    std_dev = 0.0d0
    do i = 1, count
        std_dev = std_dev + (list_plaquette(i) - mean_plaquette)**2
    end do

    deviation = sqrt((1.0d0 / (count - 1.0d0)) * std_dev)


    write(u_susc,*) beta, mean_plaquette, deviation, susceptibility              ! file with plaquette avg, std deviation and susceptibility


end subroutine metropolis_gauge








end module metropolis



