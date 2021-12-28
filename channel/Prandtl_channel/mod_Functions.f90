module Functions
implicit none
contains
!************************************************************************************************   
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
	  

!----------------------------profiles-------------------------	  
	  
	  Subroutine Profile (U, y, x, xplane)
	  
	  Real, intent(in) 	:: U(:,:), y(:,:), x(:,:), xplane
	  integer			:: i(1), j, jmax, io
	  character(len=30) :: filename, temp
	  
	  i = findloc(x(:,1), xplane)
	  
	  write(temp,'(es8.2)') x(i,1)
	  filename = 'prof_x=' // trim(adjustl(temp)) // '.plt'
	  
	  jmax = size(U,2)
	  
	  open (newunit = io, file = filename)
	  
	  Write(IO,*) 'VARIABLES = "U", "y", x' 
	  Write(IO,*) 'ZONE I=',jmax,', J=',1, ', DATAPACKING=BLOCK'
      Write(IO,'(100E25.16)') U(i,1:jmax) 
      Write(IO,'(100E25.16)') y(i,1:jmax)
	  Write(IO,'(100E25.16)') x(i,1:jmax)
	  
	  close (io)

	  end Subroutine
	  
!--------------------------------------------------	  
!----------------functions-------------------------
!----------------------------------------------------
	  
!------------------Friction coefficient--------------------

	   pure function tw (mu, u, dy)
		real, intent(in)    :: mu, dy, u(1:,1:)
		real				:: tw(size(u,1))
		integer				:: NJ
		
		NJ = size(u,2)
		tw = mu*(u(2:,NJ-2)-4*u(2:,NJ-1)+3*u(2:,NJ))*0.5/dy
	  end function
	  
!------------------Pressure calculation--------------------	  
	  pure function press (D, F, R0, U)
		integer				:: NJ
		real				:: press, G, sumF
		real, intent(in)    :: D(:), F(:), R0(:), U(:)
		press = 0.0
		G = 0.0
		NJ = size(D)
		G = sum(R0(2:NJ)*U(2:NJ) + R0(1:NJ-1)*U(1:NJ-1))
		press = sum(D(2:NJ) + D(1:NJ-1)) - G
		
		sumF = sum(F(2:NJ)+F(1:NJ-1))
		press = press/sumF
	  end function
	  
!------------------Equation of state--------------------	 
	  
	  pure function Rho (P, R0, P0, k)
		real				:: Rho
		real, intent(in)    :: P, R0, P0, k
		Rho = R0*(P/P0)**(1/k)
	  end function
		
	  
end module