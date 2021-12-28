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
	  
	  pure function tw (mu, u, dy)
		real, intent(in)    :: mu, dy, u(1:,1:)
		real				:: tw(size(u,1))
		
		tw = mu*(-3*u(2:,1)+4*u(2:,2)-u(2:,3))*0.5/dy
	  end function tw
		
	  
end module