module Navier
IMPLICIT NONE
contains
!************************************************************************************************      
       SUBROUTINE NavierStokes(P_c,U_c,V_c, dx, dy, CFL, U0, P0, NU, Uref, NITER, EPS)
	   
	   use Functions, only : boundary, strm, unstrm, alg, Solve_Tridiagonal_Matrix
	   
	   real, intent (inout) :: P_c(0:,0:),U_c(0:,0:),V_c(0:,0:)
	   real, intent (in)    :: dx, dy, CFL, U0, P0, NU, Uref, EPS
	   integer, intent (in) :: NITER
	   real, allocatable    :: Res_P(:,:),Res_U(:,:),Res_V(:,:)
	   real					:: dt, conv, visc, press
	   real, allocatable    :: Usubj(:), Vsubj(:), Psubj(:), Usubi(:), Vsubi(:), Psubi(:)
	   real					:: errorV, errorU, errorP
	   real, allocatable 	:: A(:), B(:), C(:), D(:), dV(:,:), dU(:,:), dP(:,:)
	   integer				:: t, i, j, imax, jmax, io, iu
	   
	   imax = size(P_c,1)-2
	   jmax = size(P_c,2)-2
	   
	   
	   
	   allocate (Res_P(imax,jmax),Res_U(imax,jmax),Res_V(imax,jmax)) ! internal
	   allocate (Usubj(2), Vsubj(2), Psubj(2), Usubi(2), Vsubi(2), Psubi(2)) ! internal
	   allocate (dP(imax,jmax),dU(imax,jmax),dV(imax,jmax))
	   	   
	   dt = CFL * min(0.25*dx*dx/NU, 0.25*dy*dy/NU, dx/(U0+Uref))
	   
!--------------------Set boundary condition---------------------
	   call boundary(P_c,U_c,V_c, U0, P0)

	   conv = 0.0
	   visc = 0.0
	   
	   open(newunit = io, file = 'residuals.txt')
	   do t = 1, NITER
	   
	   Res_P(:,:) = 0.0; Res_U(:,:) = 0.0; Res_V(:,:) = 0.0
			
!разбиваем область на i-столбцы, которые распределяются для параллелизма	
			do concurrent (i = 1:imax) local(Usubj, Vsubj, Psubj, conv, visc, press, j, i)
!------------------------------------------streams from previous layer--------------------
			
				Usubj(1) = strm(V_c(i,0),V_c(i,1),U_c(i,0),U_c(i,1))   !bot
				Vsubj(1) = strm(V_c(i,0),V_c(i,1),V_c(i,0),V_c(i,1))
				Psubj(1) = unstrm(V_c(i,0),V_c(i,1),P_c(i,0),P_c(i,1))
			
				do j = 1, jmax

!-----------------------------------------streams for next lay--------------------	

					Usubj(2) = strm(V_c(i,j),V_c(i,j+1),U_c(i,j),U_c(i,j+1))   !top
					Vsubj(2) = strm(V_c(i,j),V_c(i,j+1),V_c(i,j),V_c(i,j+1))
					Psubj(2) = unstrm(V_c(i,j),V_c(i,j+1),P_c(i,j),P_c(i,j+1))
															
!------------------------------------------continuity equation---------------------------
					
					if (j == 1) then
						Res_P(i,j) = -(Uref**2)*dt*(Vsubj(2)/dy) + Res_P(i,j)
					else
						Res_P(i,j) = -(Uref**2)*dt*(Vsubj(2)-Vsubj(1))/dy + Res_P(i,j)
					end if
				
						
!--------------------------------------------X Dynamics--------------------------------
										
					conv = (alg(V_c(i,j),V_c(i,j+1))*Usubj(2) - alg(V_c(i,j-1),V_c(i,j))*Usubj(1))/dy
							
					visc = NU*(U_c(i,j+1)-2*U_c(i,j)+U_c(i,j-1))/(dy**2)
					
					press = 0.0
					
					Res_U(i,j) = dt*(visc - conv - press) + Res_U(i,j)
					
!--------------------------------------------Y Dynamics--------------------------------
	
					conv = (alg(V_c(i,j),V_c(i,j+1))*Vsubj(2) - alg(V_c(i,j-1),V_c(i,j))*Vsubj(1))/dy
							
					visc = NU*(V_c(i,j+1)-2*V_c(i,j)+V_c(i,j-1))/(dy**2)
					
					press = (Psubj(2)-Psubj(1))/dy
					
					Res_V(i,j) = dt*(visc - conv - press) + Res_V(i,j)
											
!--------------------------------------------Save stream-----------------------------
					
					Usubj(1) = Usubj(2)
					Vsubj(1) = Vsubj(2)
					Psubj(1) = Psubj(2)
					
					end do

			end do
			
			
!разбиваем область на j-строки, которые распределяются для параллелизма	
			do concurrent (j = 1:jmax) local(Usubi, Vsubi, Psubi, conv, visc, press, j, i)
!------------------------------------------streams from previous layer--------------------
			
				Usubi(1) = strm(U_c(0,j),U_c(1,j),U_c(0,j),U_c(1,j))   !left
				Vsubi(1) = strm(U_c(0,j),U_c(1,j),V_c(0,j),V_c(1,j))
				Psubi(1) = unstrm(U_c(0,j),U_c(1,j),P_c(0,j),P_c(1,j))
			
				do i = 1, imax

!-----------------------------------------streams for next lay--------------------	
					
					Usubi(2) = strm(U_c(i,j),U_c(i+1,j),U_c(i,j),U_c(i+1,j))   !right
					Vsubi(2) = strm(U_c(i,j),U_c(i+1,j),V_c(i,j),V_c(i+1,j))
					Psubi(2) = unstrm(U_c(i,j),U_c(i+1,j),P_c(i,j),P_c(i+1,j))
					
!------------------------------------------continuity equation---------------------------
					
					Res_P(i,j) = -(Uref**2)*dt*(Usubi(2) - Usubi(1))/dx + Res_P(i,j)
						
!--------------------------------------------X Dynamics--------------------------------
										
					conv = (alg(U_c(i,j),U_c(i+1,j))*Usubi(2) - alg(U_c(i-1,j),U_c(i,j))*Usubi(1))/dx
							
					visc = NU*(U_c(i+1,j)-2*U_c(i,j)+U_c(i-1,j))/(dx**2)
					
					press = (Psubi(2)-Psubi(1))/dx
					
					Res_U(i,j) = dt*(visc - conv - press) + Res_U(i,j)
					
!--------------------------------------------Y Dynamics--------------------------------
	
					conv = (alg(U_c(i,j),U_c(i+1,j))*Vsubi(2) - alg(U_c(i-1,j),U_c(i,j))*Vsubi(1))/dx
							
					visc = NU*(V_c(i+1,j)-2*V_c(i,j)+V_c(i-1,j))/(dx**2)
					
					press = 0.0
					
					Res_V(i,j) = dt*(visc - conv - press) + Res_V(i,j)
											
!--------------------------------------------Save stream-----------------------------
					
					Usubi(1) = Usubi(2)
					Vsubi(1) = Vsubi(2)
					Psubi(1) = Psubi(2)
					
					end do

			end do
			
!-------------------------------------------------------------------------------------
!----------------------------------------Splitting------------------------------------
!-------------------------------------------------------------------------------------

!----------------------------------------First step-----------------------------------

			allocate (A(0:imax+1), B(0:imax+1), C(0:imax+1), D(0:imax+1))
			
			do j = 1, jmax
			
!---------------------------------------- dv* ----------------------------------------

				do i = 1, imax
					if (U_c(i,j) >= 0.0) A(i) = -dt/dx*U_c(i-1,j) - dt*NU/dx**2
					if (U_c(i,j) < 0.0) A(i) = -dt*NU/dx**2
					
					if (U_c(i,j) >= 0.0) B(i) = 1 + dt/dx*U_c(i,j) + 2*dt*NU/dx**2
					if (U_c(i,j) < 0.0) B(i) = 1 - dt/dx*U_c(i,j) + 2*dt*NU/dx**2
					
					if (U_c(i,j) >= 0.0) C(i) = -dt*NU/dx**2
					if (U_c(i,j) < 0.0) C(i) = dt/dx*U_c(i+1,j) - dt*NU/dx**2
				end do
				
				D = [0.0, (Res_V(i,j), i=1,imax), 0.0]
				
				A(0) = 0.0; A(imax+1) = 1.0
				B(0) = 1.0; B(imax+1) = -1.0
				C(0) = 0.0; C(imax+1) = 0.0
				
				call Solve_Tridiagonal_Matrix(A, B, C, D) ! solve V-dynamics equation
				dV(:,j) = D(1:imax)
				
!---------------------------------------- du* ----------------------------------------

				do i = 1, imax
					if (U_c(i,j) >= 0.0) A(i) = -(Uref*dt)**2/dx**2 - 2*dt*U_c(i-1,j)/dx - dt*NU/dx**2
					if (U_c(i,j) < 0.0) A(i) = -(Uref*dt)**2/dx**2 - dt*NU/dx**2
					
					if (U_c(i,j) >= 0.0) B(i) = 1 + 2*(Uref*dt)**2/dx**2 + 2*dt*U_c(i,j)/dx + 2*dt*NU/dx**2
					if (U_c(i,j) < 0.0) B(i) = 1 + 2*(Uref*dt)**2/dx**2 - 2*dt*U_c(i,j)/dx + 2*dt*NU/dx**2
					
					if (U_c(i,j) >= 0.0) C(i) = -(Uref*dt)**2/dx**2 - dt*NU/dx**2
					if (U_c(i,j) < 0.0) C(i) = -(Uref*dt)**2/dx**2 + 2*dt*U_c(i+1,j)/dx - dt*NU/dx**2
					
					if (U_c(i,j) >= 0.0) D(i) = Res_U(i,j)-dt/dx*(Res_P(i+1,j)-Res_P(i,j))
					if (U_c(i,j) < 0.0) D(i) = Res_U(i,j)-dt/dx*(Res_P(i,j)-Res_P(i-1,j))
				end do
				
				A(0) = 0.0; A(imax+1) = 1.0
				B(0) = 1.0; B(imax+1) = -1.0
				C(0) = 0.0; C(imax+1) = 1.0
				D(0) = 0.0; D(imax+1) = 0.0
			
				call Solve_Tridiagonal_Matrix(A, B, C, D) ! solve U-dynamics equation
				dU(:,j) = D(1:imax)
				
!---------------------------------------- dp* ----------------------------------------
				
				D(0) = 0.0; D(imax+1) = 0.0; !D(:) = dU(:,j)
				do i = 1, imax
					If (U_c(i,j) >= 0) dP(i,j) = Res_P(i,j) - dt*(Uref**2)*(D(i) - D(i-1))/dx 
					If (U_c(i,j) < 0) dP(i,j) = Res_P(i,j) - dt*(Uref**2)*(D(i+1) - D(i))/dx 
				end do
				
			end do
			deallocate (A, B, C, D)
			
!----------------------------------------Second step-----------------------------------

			allocate (A(0:jmax+1), B(0:jmax+1), C(0:jmax+1), D(0:jmax+1))
			
			do i = 1, imax
			
!---------------------------------------- du ----------------------------------------
				
				do j = 1, jmax
					if (V_c(i,j) >= 0.0) A(j) = -dt*V_c(i,j-1)/dy - dt*NU/dy**2
					if (V_c(i,j) < 0.0) A(j) = -dt*NU/dy**2
					
					if (V_c(i,j) >= 0.0) B(j) = 1 + dt*V_c(i,j)/dy + 2*dt*NU/dy**2
					if (V_c(i,j) < 0.0) B(j) = 1 - dt/dy*V_c(i,j) + 2*dt*NU/dy**2
					
					if (V_c(i,j) >= 0.0) C(j) = -dt*NU/dy**2
					if (V_c(i,j) < 0.0) C(j) = dt/dy*V_c(i,j+1) - dt*NU/dy**2
				end do
				
				D = [0.0, (dU(i,j), j=1,jmax), 0.0]
				
				A(0) = 0.0; A(jmax+1) = 1.0
				B(0) = 1.0; B(jmax+1) = -1.0
				C(0) = 1.0; C(jmax+1) = 0.0
			
				call Solve_Tridiagonal_Matrix(A, B, C, D) ! solve U-dynamics equation
				dU(i,:) = D(1:jmax)
				
!---------------------------------------- dv ----------------------------------------

				do j = 1, jmax
					if (V_c(i,j) >= 0.0) A(j) = -(Uref*dt)**2/dy**2 - 2*dt*V_c(i,j-1)/dy - dt*NU/dy**2
					if (V_c(i,j) < 0.0) A(j) = -(Uref*dt)**2/dy**2 - dt*NU/dy**2
					
					if (V_c(i,j) >= 0.0) B(j) = 1 + 2*(Uref*dt)**2/dy**2 + 2*dt*V_c(i,j)/dy + 2*dt*NU/dy**2
					if (V_c(i,j) < 0.0) B(j) = 1 + 2*(Uref*dt)**2/dy**2 - 2*dt*V_c(i,j)/dy + 2*dt*NU/dy**2
					
					if (V_c(i,j) >= 0.0) C(j) = -(Uref*dt)**2/dy**2 - dt*NU/dy**2
					if (V_c(i,j) < 0.0) C(j) = 2*dt*V_c(i,j+1)/dy - (Uref*dt)**2/dy**2 - dt*NU/dy**2
					
					if (V_c(i,j) >= 0.0) D(j) = dV(i,j) - dt/dy*(dP(i,j+1)-dP(i,j))
					if (V_c(i,j) < 0.0) D(j) = dV(i,j) - dt/dy*(dP(i,j)-dP(i,j-1))
				end do
				
				A(0) = 0.0; A(jmax+1) = 1.0
				B(0) = 1.0; B(jmax+1) = -1.0
				C(0) = 1.0; C(jmax+1) = 0.0
				D(0) = 0.0; D(jmax+1) = 0.0
			
				call Solve_Tridiagonal_Matrix(A, B, C, D) ! solve V-dynamics equation
				dV(i,:) = D(1:jmax)
				
				
!---------------------------------------- dp ----------------------------------------

				D(0) = 0.0; D(jmax+1) = 0.0; !D(:) = dV(i,:)
				do j = 1, jmax
					If (V_c(i,j) >= 0) dP(i,j) = dP(i,j) - dt*(Uref**2)*(D(j) - D(j-1))/dx 
					If (V_c(i,j) < 0) dP(i,j) = dP(i,j) - dt*(Uref**2)*(D(j+1) - D(j))/dx 
				end do
				
			end do
			
			deallocate (A, B, C, D)			

!-----------------------------------------------Residuals----------------------------------
			
			Res_U = dU
			Res_V = dV
			Res_P = dP
			
			errorU = maxval(abs(Res_U))/dt
			errorV = maxval(abs(Res_V))/dt
			errorP = maxval(abs(Res_P))/dt
			If (mod(t,100) == 0) then
				write (*,*) 'iteration: ', t
				write (*,*) 'errorU: ', errorU,' errorV: ', errorV,' errorP: ', errorP
				write (*,*) '----------------------------------------------------------'
			end if
			write (io,'(I10,3E25.16)') t, errorP, errorU, errorV
			
			If (any(isnan(U_c)).or.any(isnan(V_c)).or.any(isnan(P_c))) exit
			U_c(1:imax,1:jmax)=Res_U(1:imax,1:jmax) + U_c(1:imax,1:jmax)
			V_c(1:imax,1:jmax)=Res_V(1:imax,1:jmax) + V_c(1:imax,1:jmax)
			P_c(1:imax,1:jmax)=Res_P(1:imax,1:jmax) + P_c(1:imax,1:jmax)
			
			call boundary(P_c,U_c,V_c, U0, P0)
			If ((errorU<EPS).and.(errorV<EPS).and.(errorP<EPS)) exit
			
			
		end do
		close(io)
	   

       END  SUBROUTINE

!************************************************************************************************
end module
                    