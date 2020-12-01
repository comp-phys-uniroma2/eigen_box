program eigen
  implicit none

  integer, parameter :: dp=8
  real(dp), parameter :: Pi=3.1415926535897932846_dp
  real(dp), parameter :: eps=1e-10
  real(dp), parameter :: tol=1e-6
  real(dp), allocatable :: H(:,:), DD(:), EE(:)
  real(dp) :: dx, L, L0, E0, U0, yn, s 
  integer :: i, N, n1, n2, stencil, Nbound

  ! Variabili LAPACK 
  integer :: Nb, info, lwork, liwork
  real(dp), allocatable :: work(:)
  real(dp), allocatable :: iwork(:)
  real(dp), allocatable :: tau(:)
  integer, external :: ilaenv

  ! 1          n1       n1+n2        N
  ! ************        **************  U0
  !            *        *   
  !            *        *
  !            *        *
  !            *        *
  !            **********               0
  !            <----L---> 
  ! <------------L0------------------>

  ! Numero di punti usati per discretizzare d^2/dx^2:
  stencil = 5 

  ! hbar^2/2m in [eV nm^2]
  E0 = 0.0381_dp

  ! Lunghezza della well L e numero di suddivisioni, n2
  L = 3.0_dp
  n2 = 300
  dx = L/real(n2,dp)

  ! Altezza della well:
  U0 = 8.0_dp !eV

  ! L0
  L0 = 6.0_dp

  ! Numero di intervalli dx nelle barriere
  n1 = nint((L0-L)/2.0_dp/dx)
  print*,'n1=', n1

  ! numero totale di punti
  N = n2+2*n1
  print*,'N=',N

  ! Ridefinizione corretta di L0 come multiplo di dx
  L0=N*dx
  print*,'L0=',L0

  ! Costante hbar^2/(2 m dx**2)  [eV]
  E0 = E0/dx**2

  ! Allocazione array
  allocate(H(0:N,0:N))
  H = 0.0_dp

  ! Energia Cinetica
  select case(stencil)
  case(3)
    H(0,0) = 2.0_dp * E0
    H(0,1) = -E0
    do i = 1, N-1
      H(i,i-1) = -E0
      H(i,i) = 2.0_dp*E0 
      H(i,i+1) = -E0 
    end do    
    H(N,N-1) = -E0
    H(N,N) = 2.0_dp * E0
  case(5)
    H(0,0) =  30.0_dp/12.0_dp * E0
    H(0,1) = -16.0_dp/12.0_dp * E0
    H(0,2) =   1.0_dp/12.0_dp * E0
    H(1,0) = -16.0_dp/12.0_dp * E0
    H(1,1) =  30.0_dp/12.0_dp * E0
    H(1,2) = -16.0_dp/12.0_dp * E0
    H(1,3) =   1.0_dp/12.0_dp * E0
    do i = 2,N-2
       H(i, i-2) =   1.0_dp/12.0_dp * E0
       H(i, i-1) = -16.0_dp/12.0_dp * E0
       H(i, i  ) =  30.0_dp/12.0_dp * E0
       H(i, i+1) = -16.0_dp/12.0_dp * E0
       H(i, i+2) =   1.0_dp/12.0_dp * E0
    end do   
    H(N-1,N-3) =   1.0_dp/12.0_dp * E0 
    H(N-1,N-2) = -16.0_dp/12.0_dp * E0
    H(N-1,N-1) =  30.0_dp/12.0_dp * E0 
    H(N-1,N) =   -16.0_dp/12.0_dp * E0
    H(N,N-2) =     1.0_dp/12.0_dp * E0
    H(N,N-1) =   -16.0_dp/12.0_dp * E0
    H(N,N) =      30.0_dp/12.0_dp * E0
  case default
    stop 'ERROR: stencil 3 or 5'    
  end select


  ! Potenziale
  ! NOTA: il conto viene piu preciso se si include solo 1 dei due estremi 
  ! del bordo well. Penso perche' in questo modo si ha  <L> = L
  !
  ! *   *   *   *                       *   *   *            
  !
  !              <-------------------->
  !
  ! -   -   -   *   *   *   *   *   *   -   -   -   -   -   -
  ! |dx |   |   |                   |
  !            n1                  n1+n2
  !
  do i = 0, N
    if (i<n1 .or. i>=n1+n2) then
       H(i,i) = H(i,i) + U0
    end if   
  end do    

  ! DIAGONALIZZAZIONE 
  Nb = ilaenv(1, 'DSYTRD', 'U', N+1, -1, -1, -1)
  lwork = (N+1)*Nb
  allocate(work(lwork))
  allocate(DD(N+1))
  allocate(EE(N+1))
  allocate(tau(N+1))

  ! Tri-diagonalization
  call dsytrd('U', N+1, H, N+1, DD, EE, tau, work, lwork, info)

  deallocate(work)
 
  ! Solver QR, solo autovalori (no workspace needed)
  ! allocate(work(1))
  ! call dsteqr('N', N, DD, EE, H, N, work, info)

  ! Solver Divide&Conquer, solo autovalori (no workspace needed)
  lwork = 1
  allocate(work(lwork))
  liwork = 1
  allocate(iwork(liwork))
  call dstedc('N', N, DD, EE, H, N+1, work, lwork, iwork, liwork, info)

  E0=E0*dx**2
  ! Lista autovalori
  s = sqrt(U0*L**2/(4.0_dp*E0))

  Nbound = int(2.0_dp * s/Pi)
  print*, 'Nbound=',Nbound
 
  write(*,*) 
  write(*,'(3x,a,10x,a,13x,a,13x,a,13x,a)') '#','LAPACK','EXACT','APPROX','INF WELL' 
  write(*,'(90("-"))') 

  do i = 1, Nbound
    yn = exact(i, Pi*i/2.0_dp-eps, s, tol)
    write(*,'(i4,4(ES20.10))') i, DD(i), 4.0_dp * E0 * yn**2/L**2, &
        &  E0* Pi**2 * real(i**2,dp)/L**2/(1.0_dp+2.0_dp/L*sqrt(E0/U0))**2, &
        &  E0* Pi**2 * real(i**2,dp)/L**2
  end do

  contains

  ! ------------------------------------------------------------------------      
  ! Newton iteration to solve 
  !  x*tan(x)=sqrt(ss**2 - x**2)     n odd
  ! -x/tan(x)=sqrt(ss**2 - x**2)     n even
  function exact(n, bb, ss, tol) result(x)
    integer, intent(in) :: n    
    real(dp), intent(in) :: bb, ss
    real(dp), intent(in) :: tol
    real(dp) :: x

    real(dp) :: y, y1, error

    ! soluzioni dispari
    if (mod(n,2) == 1) then
       x=bb
       error = 1.0_dp 
       do while(error > tol)   
          y = x*tan(x)-sqrt(ss**2 - x**2) 
          y1 = tan(x)+x/(cos(x))**2+x/sqrt(ss**2 - x**2)
          x = x - y/y1
          error = abs(y)
       end do
    else
       x=bb
       error = 1.0_dp 
       do while(error > tol)   
          y = -x/tan(x)-sqrt(ss**2 - x**2) 
          y1 = -1.0_dp/tan(x) + x/(sin(x))**2+x/sqrt(ss**2 - x**2)
          x = x - y/y1
          error = abs(y)
       end do
    end if      
  
  end function exact


end program eigen

