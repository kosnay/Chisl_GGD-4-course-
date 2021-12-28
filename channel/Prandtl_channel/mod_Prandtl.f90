Module Solver
use Functions, only : Solve_Tridiagonal_Matrix, press, Rho
implicit none
contains
!************************************************************************************************      
      SUBROUTINE Prandtl(U, V, P, R, dx, dy, Mu, R0, P0, k, smax, eps)
	  
	  Real, intent(inout) :: U(:,:), V(:,:), P(:,:), R(:,:)
	  Real, intent(in)    :: dx, dy, Mu, R0, P0, k, eps
	  Integer, intent(in) :: smax
	  Real				  :: errorU, errorV, errorP, Gs, G
	  Real				  :: Cf
	  Real, allocatable   :: Un(:), Vn(:), Pn(:), A(:), B(:), C(:), D(:), F(:) !j-plate
	  Integer			  :: NI, NJ, i, j, s, io, ires, sres, iu
	  
	  NI = size(U,1); NJ = size(U,2)
	  
	  allocate(Un(NJ), Vn(NJ), Pn(NJ), A(NJ), B(NJ), C(NJ), D(NJ), F(NJ))
	  
	  Gs = 0.0
	  i=1
	  sres = 0
	  ires = 1
	  Cf = Mu/(dy*dy*R0)
	  
	  Gs = sum(R(i,2:NJ)*U(i,2:NJ)) + sum(R(i,1:NJ-1)*U(i,1:NJ-1))
	  open (newunit = io, file = 'G.txt')
	  open (newunit = iu, file = 'residuals.txt')
	  write(io,*) 'G= ', Gs
	
	  do i = 2, NI
		U(i,:) = U(i-1,:) !initial approximation from the previous layer
		V(i,:) = V(i-1,:)
		P(i,:) = P(i-1,:)
		R(i,:) = R(i-1,:)
		write (*,*) 'i plane: ', i
		write (*,*) '-----------------------------------------------'

		G = sum(R(i-1,2:NJ)*U(i-1,2:NJ)) + sum(R(i-1,1:NJ-1)*U(i-1,1:NJ-1))
		If (G /= Gs) then
			write(io,*) 'i plane: ', i
			write(io,*) 'G= ', G
			!stop
		end if
		
		s = 1
		errorU = 1.e5; errorV = 1.e5; errorP = 1.e5
		do while ((s <= smax).and.((errorU > eps).or.(errorV > eps).or.(errorP > eps)))
			A = [0.0, (-V(i,j-1)/(2*dy) - Cf, j=2,NJ-1), 0.0]
            B = [1.0, (U(i,j)/dx + 2*Cf, j=2,NJ-1), 1.0]
            C = [-1.0, (V(i,j+1)/(2*dy) - Cf, j=2,NJ-1), 0.0]
            D = [0.0, (U(i-1,j)*U(i-1,j)*R(i-1,j)/dx + P(i-1,j)/dx,j=2,NJ-1), 0.0]
			F = [0.0, (1.0/dx,j=2,NJ-1), 0.0]
			
		
			call Solve_Tridiagonal_Matrix(A, B, C, D) ! solve dynamics equation
			call Solve_Tridiagonal_Matrix(A, B, C, F) ! solve dynamics equation for F
			
			! do j = 1, NJ
				! print*, F(j), D(j)
			! end do
			!read (*,*)
			
			Pn(1) = press(D, F, R(i-1,:), U(i-1,:))
			Pn(:) = Pn(1)
			!write (*,*) 'P= ', Pn(1)
			
			R(i,1) = Rho(Pn(1), R0, P0, k) !state equation
			R(i,:) = R(i,1)
			!write(*,*) 'R= ', R(i,1)
			
			Un(:) = (D(:)-F(:)*Pn(:))/R(i,:)
			
			Vn(1) = 0.0		  				 ! solve continuity equation
			do j = 2, NJ-1                        
				Vn(j) = Vn(j-1) - dy/(2*dx*R(i,j)) * (R(i,j)*Un(j) - R(i-1,j)*U(i-1,j)+ R(i,j-1)*Un(j-1) - R(i-1,j-1)*U(i-1,j-1))
			end do
			Vn(NJ) = 0.0
			
			errorU = maxval(abs(Un(:) - U(i,:)))/maxval(abs(Un(:)))
			errorV = maxval(abs(Vn(:) - V(i,:)))/maxval(abs(Vn(:)))
			errorP = maxval(abs(Pn(:) - P(i,:)))/maxval(abs(Pn(:)))
			
			If (mod(s,100) == 0) then
				write (*,*) 'iteration: ', s
				write (*,*) 'errorU: ', errorU,' errorV: ', errorV,'errorP: ', errorP
			end if
			U(i,:) = Un(:)
			V(i,:) = Vn(:)
			P(i,:) = Pn(:)
			if (s > sres) then
				sres = s
				ires = i
			end if
			s = s + 1
			If ((s == smax).or.((errorU <= eps).and.(errorV <= eps).and.(errorP <= eps))) then
				write (iu,'(I10,3E25.16)') i, errorU, errorV, errorP
			end if
		end do
		write (*,*) '-----------------------------------------------'
	  end do
	  
	  write(*,*) 'Max iterations for',ires,' plane: ', sres

      END SUBROUTINE
	  
end module
