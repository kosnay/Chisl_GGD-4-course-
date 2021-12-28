Module Solver
use Functions, only : Solve_Tridiagonal_Matrix
implicit none
contains
!************************************************************************************************      
      SUBROUTINE Prandtl(U, V, dx, dy, Nu, smax, eps)
	  
	  Real, intent(inout) :: U(:,:), V(:,:)
	  Real, intent(in)    :: dx, dy, Nu, eps
	  Integer, intent(in) :: smax
	  Real				  :: errorU, errorV
	  Real				  :: U0, Cf
	  Real, allocatable   :: Un(:), Vn(:), A(:), B(:), C(:) !j-plate
	  Integer			  :: NI, NJ, i, j, s, sres, ires, io
	  
	  NI = size(U,1); NJ = size(U,2)
	  allocate(Un(NJ), Vn(NJ), A(NJ), B(NJ), C(NJ))
	  
	  U0 = U(1,2)
	  sres = 0
	  ires = 1
	  open (newunit = io, file = 'residuals.txt')
	  
	  do i = 2, NI
		U(i,:) = U(i-1,:) !initial approximation from the previous layer
		V(i,:) = V(i-1,:)
		write (*,*) 'i iteration: ', i
		write (*,*) '-----------------------------------------------'
		
		Cf = Nu/(dy**2)
		s = 1
		
		errorU = 1.e5; errorV = 1.e5
		do while ((s <= smax).and.((errorU > eps).or.(errorV > eps)))
		
			A = [0.0, (-V(i,j-1)/(2.0*dy) - Cf, j=2,NJ-1), 0.0]
            B = [1.0, (U(i,j)/dx + 2*Cf, j=2,NJ-1), 1.0]
            C = [0.0, (V(i,j+1)/(2*dy) - Cf, j=2,NJ-1), 0.0]
            Un = [0.0, ((U(i-1,j)**2)/dx,j=2,NJ-1), U0]
		
			call Solve_Tridiagonal_Matrix(A, B, C, Un) ! solve dynamics equation
			
			!Un(:) = D(:)
			
			Vn(1) = 0.0		  				 ! solve continuity equation
			do j = 2, NJ                         
				Vn(j) = Vn(j-1) - dy/(2*dx) * (Un(j) - U(i-1,j)+ Un(j-1) - U(i-1,j-1))
			end do
			errorU = maxval(abs(Un(:) - U(i,:)))/maxval(abs(Un(:)))
			errorV = maxval(abs(Vn(:) - V(i,:)))/maxval(abs(Vn(:)))
			write (*,*) 'iteration: ', s
			write (*,*) 'errorU: ', errorU,' errorV: ', errorV
			U(i,:) = Un(:)
			V(i,:) = Vn(:)
			if (s > sres) then
				sres = s
				ires = i
			end if
			s = s + 1
			If ((s == smax).or.((errorU <= eps).and.(errorV <= eps))) then
				write (io,'(I10,2E25.16)') i, errorU, errorV
			end if
		end do
		write (*,*) '-----------------------------------------------'
	  end do
	  
	  write(*,*) 'Max iterations for',ires,' plane: ', sres

      END SUBROUTINE
	  
end module
