module Navier
IMPLICIT NONE
contains
!************************************************************************************************      
       SUBROUTINE NavierStokes(P_c,U_c,V_c, dx, dy, CFL, U0, P0, NU, Uref, NITER, EPS)
	   
	   use Functions, only : boundary, strm, unstrm, alg, deltaUV
	   
	   real, intent (inout) :: P_c(0:,0:),U_c(0:,0:),V_c(0:,0:)
	   real, intent (in)    :: dx, dy, CFL, U0, P0, NU, Uref, EPS
	   integer, intent (in) :: NITER
	   real, allocatable    :: P_n(:,:),U_n(:,:),V_n(:,:)
	   real					:: dt, conv, visc, press
	   real, allocatable    :: Usubj(:), Vsubj(:), Psubj(:), Usubi(:), Vsubi(:), Psubi(:)
	   real					:: errorV, errorU, errorP
	   integer				:: t, i, j, imax, jmax, io
	   
		   imax = size(P_c,1)-2
		   jmax = size(P_c,2)-2
	   
	   
	   allocate (P_n(imax,jmax),U_n(imax,jmax),V_n(imax,jmax)) ! internal
	   allocate (Usubj(2), Vsubj(2), Psubj(2), Usubi(2), Vsubi(2), Psubi(2)) ! internal
	   
		write(*,*) 'Parameters:'
        write(*,*)'CFL= ', CFL, 'NI= ', imax+1, 'dx= ', dx
        write(*,*)'U0= ', U0, 'NJ= ', jmax+1, 'dy= ', dy
        write(*,*)'Uref= ', Uref, ' P0= ', P0
        !pause    
	   
	   
	   ! write (*,*) shape(P_c)
	   ! write (*,*) shape(P_n)
	   ! write (*,*) lbound(P_n)
	   ! write (*,*) ubound(P_n)
	   
	   
	   dt = CFL * min(dx, dy)/U0
	   
!--------------------Set boundary condition---------------------
	   call boundary(P_c,U_c,V_c, U0, P0)

	   conv = 0.0
	   visc = 0.0
	   
	   open(newunit = io, file = 'residuals.txt')
	   do t = 1, NITER
			
!разбиваем область на j-столбцы, которые распределяются для параллелизма	
			do i = 1,imax !local(Usubj, Vsubj, Psubj, Usubi, Vsubi, Psubi, conv, visc, press)
!------------------------------------------streams from previous layer--------------------
			
				! Usubj(1) = strm(U_c(i,0),U_c(i,1))   !bot
				! Vsubj(1) = strm(V_c(i,0),V_c(i,1))
				! Psubj(1) = unstrm(V_c(i,0),V_c(i,1),P_c(i,0),P_c(i,1))
			
				do j = 1, jmax

!-----------------------------------------streams for next lay--------------------	

					Usubj(2) = strm(U_c(i,j),U_c(i,j+1))   !top
					Vsubj(2) = strm(V_c(i,j),V_c(i,j+1))
					Psubj(2) = unstrm(V_c(i,j),V_c(i,j+1),P_c(i,j),P_c(i,j+1))
										
					Usubi(1) = strm(U_c(i-1,j),U_c(i,j))   !left
					Vsubi(1) = strm(V_c(i-1,j),V_c(i,j))
					Psubi(1) = unstrm(U_c(i-1,j),U_c(i,j),P_c(i-1,j),P_c(i,j))
					
					Usubi(2) = strm(U_c(i,j),U_c(i+1,j))   !right
					Vsubi(2) = strm(V_c(i,j),V_c(i+1,j))
					Psubi(2) = unstrm(U_c(i,j),U_c(i+1,j),P_c(i,j),P_c(i+1,j))
					
!------------------------------------------continuity equation---------------------------
					
					if (j == 1) then
						P_n(i,j) = P_c(i,j) - (Uref**2)*dt*((Usubi(2) - Usubi(1))/dx + (Vsubj(2)/dy))
					else
						P_n(i,j) = P_c(i,j) - (Uref**2)*dt*((Usubi(2) - Usubi(1))/dx + (Vsubj(2)-Vsubj(1))/dy)
					end if
					
					Usubj(1) = Usubj(2)
					Vsubj(1) = Vsubj(2)
					Psubj(1) = Psubj(2)
				end do
			end do
			errorP = maxval(abs(P_n - P_c(1:imax,1:jmax)))/(dt*t)
			P_c(1:imax,1:jmax)=P_n(1:imax,1:jmax)
			do i = 1,imax !local(Usubj, Vsubj, Psubj, Usubi, Vsubi, Psubi, conv, visc, press)
!------------------------------------------streams from previous layer--------------------
			! Usubj(1) = strm(U_c(i,0),U_c(i,1))   !bot
            ! Vsubj(1) = strm(V_c(i,0),V_c(i,1))
            ! Psubj(1) = unstrm(V_c(i,0),V_c(i,1),P_c(i,0),P_c(i,1))

			
				do j = 1, jmax

!-----------------------------------------streams for next lay--------------------
					Usubj(1) = strm(U_c(i,j-1),U_c(i,j))   !bot
					Vsubj(1) = strm(V_c(i,j-1),V_c(i,j))
					Psubj(1) = unstrm(V_c(i,j-1),V_c(i,j),P_c(i,j-1),P_c(i,j))
                  
	              	Usubj(2) = strm(U_c(i,j),U_c(i,j+1))   !top
					Vsubj(2) = strm(V_c(i,j),V_c(i,j+1))
					Psubj(2) = unstrm(V_c(i,j),V_c(i,j+1),P_c(i,j),P_c(i,j+1))
										
					Usubi(1) = strm(U_c(i-1,j),U_c(i,j))   !left
					Vsubi(1) = strm(V_c(i-1,j),V_c(i,j))
					Psubi(1) = unstrm(U_c(i-1,j),U_c(i,j),P_c(i-1,j),P_c(i,j))
					
					Usubi(2) = strm(U_c(i,j),U_c(i+1,j))   !right
					Vsubi(2) = strm(V_c(i,j),V_c(i+1,j))
					Psubi(2) = unstrm(U_c(i,j),U_c(i+1,j),P_c(i,j),P_c(i+1,j))
						
!--------------------------------------------X Dynamics--------------------------------
										
					conv = (alg(U_c(i,j),U_c(i+1,j))*Usubi(2) - alg(U_c(i-1,j),U_c(i,j))*Usubi(1))/dx &
							+ (alg(V_c(i,j),V_c(i,j+1))*Usubj(2) - alg(V_c(i,j-1),V_c(i,j))*Usubj(1))/dy
					!conv = deltaUV(U_c(i-1:i+1,j), U_c(i-1:i+1,j), dx) + deltaUV(V_c(i,j-1:j+1), U_c(i,j-1:j+1), dy) 
							
					visc = NU*((U_c(i+1,j)-2*U_c(i,j)+U_c(i-1,j))/(dx**2) + (U_c(i,j+1)-2*U_c(i,j)+U_c(i,j-1))/(dy**2))
					
					press = (Psubi(2)-Psubi(1))/dx
					
					U_n(i,j) = U_c(i,j) + dt*(visc - conv - press)
					
!--------------------------------------------Y Dynamics--------------------------------
	
					conv = (alg(U_c(i,j),U_c(i+1,j))*Vsubi(2) - alg(U_c(i-1,j),U_c(i,j))*Vsubi(1))/dx &
							+ (alg(V_c(i,j),V_c(i,j+1))*Vsubj(2) - alg(V_c(i,j-1),V_c(i,j))*Vsubj(1))/dy
					!conv = deltaUV(U_c(i-1:i+1,j), V_c(i-1:i+1,j), dx) + deltaUV(V_c(i,j-1:j+1), V_c(i,j-1:j+1), dy)
							
					visc = NU*((V_c(i+1,j)-2*V_c(i,j)+V_c(i-1,j))/(dx**2) + (V_c(i,j+1)-2*V_c(i,j)+V_c(i,j-1))/(dy**2))
					
					press = (Psubj(2)-Psubj(1))/dy
					
					V_n(i,j) = V_c(i,j) + dt*(visc - conv - press)
											
!--------------------------------------------Save stream-----------------------------
					
					! Usubj(1) = Usubj(2)
					! Vsubj(1) = Vsubj(2)
					! Psubj(1) = Psubj(2)
					
					end do

			end do
!-----------------------------------------------Residuals----------------------------------
			errorU = maxval(abs(U_n - U_c(1:imax,1:jmax)))/(dt*t)
			errorV = maxval(abs(V_n - V_c(1:imax,1:jmax)))/(dt*t)
			
			If (mod(t,1) == 0) then
				write (*,*) 'iteration: ', t
				write (*,*) 'errorU: ', errorU,' errorV: ', errorV,' errorP: ', errorP
				write (*,*) '----------------------------------------------------------'
				! write (*,*) 'U: ', U_n(1,:),' V: ', V_n(1,:),' P: ', P_n(1,:)
				! write (*,*) '----------------------------------------------------------'
			end if
			write (io,'(I10,3E25.16)') t, errorP, errorU, errorV
			
			If (any(isnan(U_c)).or.any(isnan(V_c)).or.any(isnan(P_c))) exit
			U_c(1:imax,1:jmax)=U_n(1:imax,1:jmax)
			V_c(1:imax,1:jmax)=V_n(1:imax,1:jmax)
			
			
			If ((errorU<EPS).and.(errorV<EPS).and.(errorP<EPS)) exit
			call boundary(P_c,U_c,V_c, U0, P0)
			
		end do
		close(io)
	   

       END  SUBROUTINE

!************************************************************************************************
end module
                    