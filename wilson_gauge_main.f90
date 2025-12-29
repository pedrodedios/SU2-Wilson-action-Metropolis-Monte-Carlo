!====================================================================================================
! Name        : wilson_gauge_main.f90
! Author      : Pedro de Dios
! Description :  Monte Carlo simulation of SU(2) pure gauge theory on a 4D Euclidean 
!               lattice using the Metropolis algorithm.
!   The code implements:
!     - Periodic boundary conditions on a 4D hypercubic lattice
!     - SU(2) link variables stored as 2×2 complex matrices
!     - Random SU(2) hot-start initialization (Haar measure)
!     - Local Metropolis updates using small SU(2) rotations (**optional**)
!     - Computation of staples and plaquette observables
!     - Measurement of the average plaquette and specific heat 
!     - OpenMP parallelization for plaquette measurements
!
!   The simulation scans multiple values of the inverse coupling β and
!   measures gauge-invariant observables after thermalization.
!
!   Intended for educational and research use in lattice gauge theory.
!
! Note:   !!! Edit **omp_set_num_threads(4)** according to your # of CPU's for better performance !!!
!====================================================================================================






program wilson_gauge_main




use gauge_utils          ! math and SU(2) utilities
use metropolis           ! metropolis algorithm 
use omp_lib              ! Enables OpenMP shared-memory parallelism (thread management and runtime routines)




implicit none

!call random_seed()

! Parameters
integer :: XS, YS, ZS,TS, num_sweeps, thermalization, measuring_gap
real(8) :: beta, beta_i, beta_f,t_i,t_f,t_start, t_end,eps, beta_step
integer :: n, measuring_gag, num_steps

character(len=100) :: fname_action, fname_plaq, fname_susc
integer :: u_action, u_plaq, u_susc




 ! --- Set free parameters ---
XS = 6   !X dimension
YS = 6   !Y dimension
ZS = 6   !Z dimension
TS = 4   !T dimension

num_sweeps = 10000          !# of sweeps (~80k-100k for reliable statistics)
thermalization = 3000   !# of sweeps before measuring observables (aprox 40k )


beta_i = 1.0d0            ! Initial beta
beta_f = 4.0d0            ! Final beta
beta_step = 0.2d0        ! Increase in beta

num_steps = int ((beta_f - beta_i) / beta_step)

eps = 0.2d0   ! tune this parameter to explore small rotated proposals U' = R*U  (**optional**)
    

call cpu_time(t_start)
call omp_set_num_threads(4)   !!!! EDIT ACCORDING TO # OF CPU's available!!!

open(unit=u_susc, file='susceptibility_avg', status='replace')



    do n = 0, num_steps  !Run MC-Metropolis for different couplings 
  
        beta = beta_i + beta_step*dble(n)

        print *, 'running MC-metropoplis, beta:', beta 

        call cpu_time(t_i)

        ! File names for saving data
        write(fname_action, '(A,F6.3,A)') 'action_evolution_beta_', beta, '.txt'
        write(fname_plaq,   '(A,F6.3,A)') 'plaquette_evolution_beta_', beta, '.txt'


        ! Open files with automatically assigned unit numbers
        open(newunit=u_action, file=trim(fname_action), status='replace', action='write')
        open(newunit=u_plaq,   file=trim(fname_plaq),   status='replace', action='write')



 
        ! tune measuring gab when measuring observables
        if (beta < 1.5d0) measuring_gap = 10
        if (beta >= 1.5d0 .and. beta < 1.9d0) measuring_gap = 20
        if (beta >= 1.9d0) measuring_gap = 25



        call metropolis_gauge(beta, XS, YS, ZS, TS, eps, num_sweeps, thermalization, measuring_gap, &
             u_action, u_plaq, u_susc)  ! MC-Metropolis algorithm 
  

        call cpu_time(t_f)
 
        print *, 'end beta. Running time (s)', t_f - t_i 
        print *, ' '


    end do ! beta loop



call cpu_time(t_end)
print *, "Total time (s) =", t_end - t_start







end program wilson_gauge_main







