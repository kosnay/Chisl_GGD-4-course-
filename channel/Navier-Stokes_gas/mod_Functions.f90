module Functions
implicit none
contains
!************************************************************************************************
	subroutine boundary(P_c,U_c,V_c,R_c, U0, P0, R0, k)
	
	real, intent (inout) :: P_c(0:,0:),U_c(0:,0:),V_c(0:,0:),R_c(0:,0:)
	real, intent (in)    :: U0, P0, R0, k
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
	
	!wall
	  P_c(1:imax,NJ) = P_c(1:imax,jmax)
	  V_c(1:imax,NJ) = -V_c(1:imax,jmax)
	  U_c(1:imax,NJ) = -U_c(1:imax,jmax)
	  R_c(1:imax,NJ) = R_c(1:imax,jmax)
	!sym
	  P_c(1:imax,0) = P_c(1:imax,1)
	  V_c(1:imax,0) = -V_c(1:imax,1)
	  U_c(1:imax,0) = U_c(1:imax,1)
	  R_c(1:imax,0) = R_c(1:imax,1)
	!inlet
	  P_c(0,0:NJ) = P_c(1,0:NJ)
	  V_c(0,0:NJ) = 0.0
	  U_c(0,0:NJ) = U0
	  R_c(0,0:NJ) = R_c(1,0:NJ)
	!outlet
	  P_c(NI,0:NJ) = P0
	  V_c(NI,0:NJ) = V_c(imax,0:NJ)
	  U_c(NI,0:NJ) = U_c(imax,0:NJ)
	  R_c(NI,0:NJ) = Rho(P_c(NI,0:NJ), R0, P0, k)
	   
	   	   
	end subroutine
	
!----------------------------profiles-------------------------	  
	  
	  Subroutine Profile (U, y, x, xplane)
	  
	  Real, intent(in) 	:: U(0:,0:), y(:,:), x(:,:), xplane
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
		real				:: tw(size(u,1))
		integer				:: NJ
		
		NJ = size(u,2) - 1
		tw = mu*(u(:,NJ)-u(:,NJ-1))/dy
	end function
	
	
	real pure elemental function Press(R, R0, P0, k) result(P)
		real, intent(in)    :: R, R0, P0, k
		P = P0*(R/R0)**k
	end function
	
	real pure elemental function Rho(P, R0, P0, k) result(R)
		real, intent(in)    :: P, R0, P0, k
		R = R0*(P/P0)**(1/k)
	end function
	
!************************************************************************************************
end module


