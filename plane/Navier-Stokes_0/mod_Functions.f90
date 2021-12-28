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
	
	
	pure function tw (mu, u, dy, NI)
		real			    :: tw(size(u,1))
		real, intent(in)    :: mu, dy, u(1:,0:)
		tw = mu*(u(:,1)-u(:,0))/dy
	end function
	
	
	! real pure function deltaUV (ak, a, d) ! a(1) - Аi-1,j, a(2) - Ai,j, a(3) - Ai+1,j
		! real, intent (in) :: ak(1:), a(1:), d   ! d - dx/dy
		! deltaUV = (alg(ak(2),ak(3))*strm(a(2),a(3)) - alg(ak(1),ak(2))*strm(a(1),a(2)))/d
	! end function
	
!************************************************************************************************
end module


