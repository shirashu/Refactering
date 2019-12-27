! Time-stamp: <Fri Apr 27 08:00:59 JST 2018> 
!
! 1D DRIFT-DIFFUSION SIMULATOR / Nobuya Mori (Osaka Univ)
!

!
! Constants, parameters, and conditions
!

MODULE mod_cst

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(1.0d0) ! double precision

  ! physical constants

  REAL(DP), PARAMETER :: PV = 8.854187817620d-12 ! electric const  (F/m)
  REAL(DP), PARAMETER :: EC = 1.6021766208d-19   ! elementary charge (C)
  REAL(DP), PARAMETER :: KB = 1.38064852d-23     ! Boltzmann const (J/K)

  ! material parameters and conditions

  REAL(DP), PARAMETER :: TL   = 300d0        ! lattice temperature   (K)
  REAL(DP), PARAMETER :: VT   = KB * TL / EC ! thermal voltage       (V)
  REAL(DP), PARAMETER :: NI   = 1.08d16      ! intrinsic density   (/m3)
  REAL(DP), PARAMETER :: DC   = 11.9d0 * PV  ! dielectric constant (F/m)
  REAL(DP), PARAMETER :: MUP  = 0.050d0      ! hole mobility     (m2/Vs)
  REAL(DP), PARAMETER :: MUN  = 0.147d0      ! electron mobility (m2/Vs)
  REAL(DP), PARAMETER :: TAUP = 10.0d-9      ! hole life time        (s)
  REAL(DP), PARAMETER :: TAUN = 10.0d-9      ! electron life time    (s)

END MODULE mod_cst

!
! Bernoulli function
!

MODULE mod_bernoulli

  USE mod_cst

  IMPLICIT NONE

CONTAINS

  REAL(DP) FUNCTION Bernoulli(x)

    REAL(DP), INTENT(in) :: x

    IF (ABS(x) < 1d-2) THEN
       Bernoulli = 1d0 - x * (0.5d0 - x * (1d0 / 12 - x**2 / 720))
    ELSE
       Bernoulli = x / (EXP(x) - 1d0)
    END IF

  END FUNCTION Bernoulli

END MODULE mod_bernoulli

!
! Generation and recombination rate
!

MODULE mod_gr

  USE mod_cst

  IMPLICIT NONE

CONTAINS

  REAL(DP) FUNCTION gr_rate(p, n)

    REAL(DP), INTENT(in) :: p, n ! normalized by the intrinsic density

    gr_rate = -NI * (p * n - 1d0) / (TAUP * (p + 1d0) + TAUN * (n + 1d0))

  END FUNCTION gr_rate

END MODULE mod_gr

!
! Solve tridiagonal linear systems A x = b
!   CALL tridg(a)
!     a(n,4) double precision matrix
!     a(2:n,  1) [in]  lower sub-diagonal part of A
!     a(1:n,  2) [in]  diagonal part of A
!     a(1:n-1,3) [in]  upper sub-diagonal of A
!     a(1:n,  4) [in]  right hand side, b
!                [out] solution, x
! See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
!

MODULE mod_tridg

  USE mod_cst

  IMPLICIT NONE

CONTAINS

  SUBROUTINE tridg(a)

    REAL(DP), INTENT(inout) :: a(:,:) ! a(n,4)

    INTEGER :: n, i

    IF (SIZE(a, 2) /= 4) STOP
    n = SIZE(a, 1)

    DO i = 2, n
       a(i,1) = a(i,1) / a(i-1,2)
       a(i,2) = a(i,2) - a(i,1) * a(i-1,3)
       a(i,4) = a(i,4) - a(i,1) * a(i-1,4)
    END DO

    a(n,4) = a(n,4) / a(n,2)
    DO i = n - 1, 1, -1
       a(i,4) = (a(i,4) - a(i,3) * a(i+1,4)) / a(i,2)
    END DO

  END SUBROUTINE tridg

END MODULE mod_tridg

!
! Solve Poisson equation
!

MODULE mod_poisson

  USE mod_cst
  USE mod_tridg

  IMPLICIT NONE


CONTAINS

  SUBROUTINE poisson(eta, p, n, c, psi, dpsi)
    
    REAL(DP), INTENT(in)  :: eta
    REAL(DP), INTENT(in)  :: p(0:), n(0:), c(0:), psi(0:)
    REAL(DP), INTENT(out) :: dpsi(1:)

    REAL(DP), ALLOCATABLE :: a(:,:)
    INTEGER :: m, i

    m = UBOUND(psi, 1)
    ALLOCATE(a(m-1,4))

    a = 0d0
    a(2:m-1,1) = -1d0
    a(1:m-1,2) = 2d0 + eta * (p(1:m-1) + n(1:m-1))
    a(1:m-2,3) = -1d0
    a(1:m-1,4) = eta * (p(1:m-1) - n(1:m-1) + c(1:m-1))

    DO i = 1, m - 1
       a(i,4) = a(i,4) + psi(i-1) - 2d0 * psi(i) + psi(i+1)
    END DO

    CALL tridg(a)

    dpsi(1:m-1) = a(1:m-1,4)

  END SUBROUTINE poisson

END MODULE mod_poisson

!
! Solve continuity equation
!

MODULE mod_continuity

  USE mod_cst
  USE mod_bernoulli
  USE mod_tridg

  IMPLICIT NONE


CONTAINS

  SUBROUTINE continuity(psi, gr, pn)

    REAL(DP), INTENT(in)    :: psi(0:)
    REAL(DP), INTENT(in)    :: gr(1:)
    REAL(DP), INTENT(inout) :: pn(0:)
    
    REAL(DP), ALLOCATABLE :: a(:,:)
    INTEGER :: m, i

    m = UBOUND(psi, 1)
    ALLOCATE(a(m-1,4))

    a = 0d0

    DO i = 2, m - 1
       a(i,1) = -Bernoulli(psi(i) - psi(i-1))
    END DO

    DO i = 1, m - 1
       a(i,2) = Bernoulli(psi(i+1) - psi(i)) &
              + Bernoulli(psi(i-1) - psi(i))
    END DO

    DO i = 1, m - 2
       a(i,3) = -Bernoulli(psi(i) - psi(i+1))
    END DO

    a(1:m-1,4) = gr(1:m-1)

    ! boundary condition

    a(1,4)   = a(1,4)   + Bernoulli(psi(1)   - psi(0)) * pn(0)
    a(m-1,4) = a(m-1,4) + Bernoulli(psi(m-1) - psi(m)) * pn(m)

    CALL tridg(a)

    pn(1:m-1) = a(1:m-1,4)

  END SUBROUTINE continuity

END MODULE mod_continuity

!
! Evaluate current
!

MODULE mod_current

  USE mod_cst
  USE mod_bernoulli

  IMPLICIT NONE

CONTAINS

  SUBROUTINE current(psi, pn, Jpn)
    
    REAL(DP), INTENT(in)  :: psi(0:) ! normalized by VT
    REAL(DP), INTENT(in)  :: pn(0:)
    REAL(DP), INTENT(out) :: Jpn(0:)

    INTEGER :: m, i

    m = UBOUND(psi, 1)

    DO i = 0, m - 1
       Jpn(i) = Bernoulli(psi(i) - psi(i+1)) * pn(i+1) &
              - Bernoulli(psi(i+1) - psi(i)) * pn(i)
    END DO

  END SUBROUTINE current

END MODULE mod_current

!
! Main routine
!

PROGRAM drift_diffusion

  USE mod_cst
  USE mod_gr
  USE mod_poisson
  USE mod_continuity
  USE mod_current

  IMPLICIT NONE

  REAL(DP) :: tolerance            ! potential tolerance (V)

  REAL(DP) :: Vstp                 ! applied bias step (V)
  INTEGER  :: Vnum                 ! number of bias points
  REAL(DP) :: L                    ! device length (m)
  INTEGER  :: m                    ! number of mesh points
  REAL(DP) :: Nd                   ! doping density in n-region (/m3)
  REAL(DP) :: Na                   ! doping density in p-region (/m3)

  REAL(DP), ALLOCATABLE :: p(:)    ! hole density (/m3)
  REAL(DP), ALLOCATABLE :: n(:)    ! electron density (/m3)
  REAL(DP), ALLOCATABLE :: c(:)    ! net impurity density Nd - Na (/m3)
  REAL(DP), ALLOCATABLE :: psi(:)  ! potential (V)
  REAL(DP), ALLOCATABLE :: dpsi(:) ! potential update (V)
  REAL(DP), ALLOCATABLE :: Jp(:)   ! hole current density (A/m2)
  REAL(DP), ALLOCATABLE :: Jn(:)   ! electron current density (A/m2)
  REAL(DP), ALLOCATABLE :: gr(:)   ! net GR rate (/s)
  REAL(DP), ALLOCATABLE :: gp(:)   ! normalized GR rate, Gamma_p
  REAL(DP), ALLOCATABLE :: gn(:)   ! normalized GR rate, Gamma_n

  REAL(DP) :: h, Vapp, psi0, eta, nni, ujp, ujn, pfp, pfn, residue
  INTEGER  :: loop, i, iv

  ! read conditions

  READ(*,*) L
  READ(*,*) m
  READ(*,*) Nd
  READ(*,*) Na
  READ(*,*) Vstp
  READ(*,*) Vnum
  READ(*,*) tolerance

  ! set up the system

  h = L / m ! mesh width (m)

  ALLOCATE(p(0:m), n(0:m), c(0:m), psi(0:m), dpsi(1:m-1))
  ALLOCATE(Jp(0:m-1), Jn(0:m-1), gr(1:m-1), gp(1:m-1), gn(1:m-1))

  c(0:(m-1)/2) = -Na ! p-region
  c((m+1)/2:m) =  Nd ! n-region

  ! initial guess for p & n

  p = 0.5d0 * (SQRT(c**2 + 4 * NI**2) - c)
  n = 0.5d0 * (SQRT(c**2 + 4 * NI**2) + c)

  ! constants for normalization

  eta = EC / h / DC / VT     ! for Poisson

  ujp = EC * MUP * VT / h**4 ! unit of hole current
  ujn = EC * MUN * VT / h**4 ! unit of electron current

  pfp = h**5 / (MUP * VT)    ! pre-factor for Gamma_p
  pfn = h**5 / (MUN * VT)    ! pre-factor for Gamma_n

  nni = NI * h**3            ! normalized intrinsic density

  ! normalization

  p = p * h**3
  n = n * h**3
  c = c * h**3

  ! initial guess for normalized potential

  psi = LOG(n / nni)
  psi0 = psi(0)

  ! main loop

  DO iv = 0, Vnum

     Vapp = iv * Vstp
     psi(0) = psi0 + Vapp / VT ! normalized applied voltage

     ! self-consistent loop

     loop = 0
     DO WHILE (.TRUE.)

        loop = loop + 1

        ! generation-recombination rate

        DO i = 1, m - 1
           gr(i) = gr_rate(p(i) / nni, n(i) / nni)
        END DO
        gp = pfp * gr
        gn = pfn * gr

        ! solve continuity equations

        CALL continuity( psi, gp, p)
        CALL continuity(-psi, gn, n)

        ! solve Poisson equation

        CALL poisson(eta, p, n, c, psi, dpsi)
        psi(1:m-1) = psi(1:m-1) + dpsi(1:m-1)

        ! converge?

        residue = SQRT(SUM(dpsi(1:m-1)**2) / (m - 1)) * VT ! (V)
        IF (residue < tolerance) EXIT

     END DO

     ! evaluate current

     CALL current( psi, p, Jp); Jp = -ujp * Jp
     CALL current(-psi, n, Jn); Jn =  ujn * Jn

     WRITE(*,'(4(es15.8,1x),i3)') &
          Vapp, Jp(0), Jn(0), Jp(0) + Jn(0), loop

  END DO

  ! restore the units

  p = p / h**3
  n = n / h**3
  c = c / h**3

  psi = psi * VT

  ! save results

  OPEN(110,file='dd_out_m.txt')
  DO i = 0, m
     WRITE(110,'(5(es15.8,1x))') &
          i * h, psi(i), p(i), n(i), p(i) - n(i) + c(i)
  END DO
  CLOSE(110)

  OPEN(120,file='dd_out_a.txt')
  DO i = 0, m - 1
     WRITE(120,'(5(es15.8,1x))') (i + 0.5d0) * h, &
          -(psi(i+1) - psi(i)) / h, Jp(i), Jn(i), Jp(i) + Jn(i)
  END DO
  CLOSE(120)

END PROGRAM drift_diffusion
