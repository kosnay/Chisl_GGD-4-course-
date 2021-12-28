        Program Pr
		use Solver, only : Prandtl
		use InOutput, only : Output_Fields, Output_wall
		use Functions, only : tw, Profile
        Implicit none
		
		
		character(len=100),parameter :: input_file='input_params.nml'
        INTEGER 		 			 :: NI, NJ, NITER, SMAX
        INTEGER 		 			 :: I, J, IO
        REAL 			 			 :: L,H,U0,MU,Nu,R0,P0,k
        REAL 			 			 :: dx,dy,CFL,EPS
        REAL,ALLOCATABLE 			 :: X_Node(:,:),Y_Node(:,:), Cf(:), Cft(:)
        REAL,ALLOCATABLE 			 :: U_n(:,:),V_n(:,:),P_n(:,:),R_n(:,:)
		
		namelist /input_params/ NI,NJ,L,H,U0,MU,R0,P0,k,EPS,SMAX
		open(newunit = IO, file = input_file, action = 'read')
		read(IO, nml = input_params)
		close(IO)
		
		write(*,*) 'Read input file - Done'
		dx=L/(NI-1)
        dy=H/(NJ-1)
		
		NU=MU/R0
   
        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
		

!----------------- Node variables -----------------------------
        allocate(U_n(NI,NJ))  ! Velocity U
        allocate(V_n(NI,NJ))  ! Velocity V
        allocate(P_n(NI,NJ))  ! Pressure
		allocate(R_n(NI,NJ))  ! Density
		allocate(Cf(2:NI),Cft(2:NI))    ! Coefficent of friction

!----------------- Coordinate of nodes ------------------------
        

        DO I=1,NI
          DO J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          END DO
        END DO

!----------------- Parameters ------------------------

        

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
		write(*,*)'k= ', k
        write(*,*)'ReL= ', U0*2*H/NU
		read (*,*) 

!----------------- Initial fields -----------------------------

        DO I=1,NI
          DO J=1,NJ
            U_n(I,J)=U0
            V_n(I,J)=0.0
            P_n(I,J)=P0
			R_n(I,J)=R0
          ENDDO
        ENDDO

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Prandtl equations'      
        call  Prandtl(U_n, V_n, P_n, R_n, dx, dy, Mu, R0, P0, k, smax, eps)
		
		Cf(2:) = -tw(mu, U_n(:,:), dy)*2/(R0*U0*U0) ! "-" - normal to wall
		Cft(2:) = 12.0/(U0*2*H/NU)

 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data' 
        Open(newunit = IO,FILE='data_pr.plt')
        Call Output_Fields(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
        Close(IO)
		write(*,*) 'Output wall data' 
		call Output_wall(X_Node(2:,1),Cf(:), Cft(:))
		
!------------------profiles------------------------------------		
		write(*,*) 'Output profiles data' 
		j = (NI-1)/10
		
		do i = 1,10,2
			write(*,*)'xplane ', X_Node(j*i+1,1)
			call Profile(U_n, Y_Node, X_Node, X_Node(j*i+1,1))
		end do
     
        END PROGRAM
   
