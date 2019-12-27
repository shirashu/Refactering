! Monte Carlo Simulation 
!
!
!

!
! Constants, parameters, and conditions
!

module mod_cst
    implicit none

    integer, parameter :: DP = KIND(1.0d0)  ! double precision

    ! physical contants 

    real(DP), parameter :: EP0 = 8.854187817620d-12 ! electric const  (F/m)
    real(DP), parameter :: Q = 1.6021766208d-19     ! elementary charge (C)
    real(DP), parameter :: KB = 1.38064852d-23      ! Boltzmann const (J/K)
    real(DP), parameter :: h = 1.05459d-34          ! dhirach const (Js)　=hbar
    real(DP), parameter :: am0 = 9.10953d-31        ! electron math (kg)
    real(DP), parameter :: PI = 3.1415926           ! pi

    ! material parameters and conditions

    REAL(DP), PARAMETER :: TL   = 300d0             ! lattice temperature   (K)


contains
    
end module mod_cst