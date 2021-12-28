module Functions
implicit none
contains
!************************************************************************************************
	subroutine boundary(P_c,U_c,V_c, U0, P0)
	
	real, intent (inout) :: P_c(0:,0:),U_c(0:,0:),V_c(0:,0:)
	real, intent (in)    :: U0, P0
	integer				 :: imax, jmax, NI, NJ
	
	imax = size(P_c,1)-2
	jmax = size(P_c,2)-2
	
	NI = imax + 1
	NJ = jmax + 1
	
	! write (*,*) strm (1.0,-2.0)
	! stop
	
	! write (*,*) shape(P_c)
	! write (*,*) lbound(P_c)
	! write (*,*) ubound(P_c)
	! write (*,*) NI, NJ
	! write (*,*) imax, jmax
	! print*,'--------------------------'
	!stop
	
	!inlet
	  P_c(0,1:jmax) = P_c(1,1:jmax)
	  V_c(0,1:jmax) = 0.0
	  U_c(0,1:jmax) = U0
	!outlet
	  P_c(NI,1:jmax) = 0.0
	  V_c(NI,1:jmax) = V_c(imax,1:jmax)
	  U_c(NI,1:jmax) = U_c(imax,1:jmax)
	!wall
	  P_c(1:imax,0) = P_c(1:imax,1)
	  V_c(1:imax,0) = -V_c(1:imax,1)
	  U_c(1:imax,0) = -U_c(1:imax,1)
	!external
	  V_c(1:imax,NJ) = V_c(1:imax,jmax)
	  where (V_c(1:imax,NJ) >= 0.0)
		P_c(1:imax,NJ) = P0
		U_c(1:imax,NJ) = U_c(1:imax,jmax)
	  elsewhere (V_c(1:imax,NJ) < 0)
		P_c(1:imax,NJ) = P_c(1:imax,jmax)
		U_c(1:imax,NJ) = U0
	  end where
	   
	   	   
	end subroutine
	
	real pure function strm (a, b, c, d) ! a - Аi,j, b - Ai+1,j
		real, intent (in) :: a, b, c, d
		if (alg(a, b) >= 0.0 ) strm = c
		if (alg(a, b) < 0.0 ) strm = d
	
	end function strm
	
	real pure function unstrm (a, b, c, d) ! a,c - Аi,j, b,d - Ai+1,j
		real, intent (in) :: a, b, c, d
		if (alg(a, b) >= 0.0 ) unstrm = d
		if (alg(a, b) < 0.0 ) unstrm = c
	
	end function unstrm
	
	real pure function alg (a, b) ! a - Аi,j, b - Ai+1,j foк boundary
		real, intent (in) :: a, b
		alg = (a + b)*0.5
	
	end function alg
	
	
	pure function tw (mu, u, dy)
		real, intent(in)	:: mu, dy, u(1:,0:)
		real			 	:: tw(size(u,1))
		
		tw = mu*(u(:,1)-u(:,0))/dy
	end function
	
	
	SUBROUTINE Solve_Tridiagonal_Matrix(A, B, C, D)
		
		Real, intent(inout) :: D(:)
		Real, intent(in)    :: A(:), B(:), C(:)
		Real, allocatable   :: alph(:), beta(:)
		Real				  :: Cf, z
		Integer			  :: i, NI
		
		NI = size(A)
		allocate (alph(NI), beta(NI))
	  
		
!--------------- прямой ход---------------
		
		alph(1) = -C(1)/B(1) !-C/B
		beta(1) = D(1)/B(1) ! D/B
		
!расчет коэф альфа
		do i=2, NI-1
			z = B(i) + A(i)*alph(i-1) !gamma coefficient
			
			alph(i) = -C(i)/z
			beta(i) = (D(i) - A(i)*beta(i-1))/z
		end do
		
		z = B(NI) + A(NI)*alph(NI-1)
		D(NI) = (D(NI) - A(NI)*beta(NI-1))/z 
	
!--------------- обратный ход ----------------
		do i=NI-1,1,-1 
			D(i)=D(i+1)*alph(i)+beta(i)
		end do
	END SUBROUTINE
	
	Subroutine Blazius (U, y, x, xplane)
	  
	  Real, intent(in) 	:: U(:,:), y(:,:), x(:,:), xplane
	  Real, allocatable :: phi(:), eta(:)
	  Real				:: mindU, maxU
	  integer			:: i(1), j, maxj, jmax, io
	  character(len=30) :: filename, temp
	  
	  i = findloc(x(:,1), xplane)
	  
	  write(temp,'(es8.2)') x(i,1)
	  filename = 'prof_x=' // trim(adjustl(temp)) // '.plt'
	  
	  maxU = maxval(U(i,:))
	  
	  jmax = size(U,2)
	  
	  do j = 1, jmax
		mindU = abs(U(i(1),j) - maxU)/maxU
		if (mindU <= 0.01) then
			maxj = j
			exit
		end if
	  end do
		
	  allocate (phi(1:maxj), eta(1:maxj))
	  phi = U(i(1),1:maxj)/U(i(1),maxj)
	  eta = y(i(1),1:maxj)/y(i(1),maxj)
	  
	  open (newunit = io, file = filename)
	  
	  Write(IO,*) 'VARIABLES = "phi", "eta", x' 
	  Write(IO,*) 'ZONE I=',maxj,', J=',1, ', DATAPACKING=BLOCK'
      Write(IO,'(100E25.16)') phi(1:maxj) 
      Write(IO,'(100E25.16)') eta(1:maxj)
	  Write(IO,'(100E25.16)') x(i,1:maxj)
	  
	  close (io)

	  end Subroutine
	
!************************************************************************************************
end module


