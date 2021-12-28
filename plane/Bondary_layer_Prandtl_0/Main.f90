        Program Pr
		use Solver, only : Prandtl
		use InOutput, only : Output_Fields, Output_wall
		use Functions, only : tw, Blazius
        Implicit none
		
		
		character(len=100),parameter :: input_file='input_params.nml'
        INTEGER 		 			 :: NI, NJ, NITER, SMAX
        INTEGER 		 			 :: I, J, IO
        REAL 			 			 :: L,H,U0,MU,Nu,R0,P0
        REAL 			 			 :: dx,dy,CFL,EPS
        REAL,ALLOCATABLE 			 :: X_Node(:,:),Y_Node(:,:), Cf(:), Cft(:)
        REAL,ALLOCATABLE 			 :: U_n(:,:),V_n(:,:),P_n(:,:),R_n(:,:)
		
		namelist /input_params/ NI,NJ,L,H,U0,MU,R0,P0,EPS,SMAX
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
        write(*,*)'ReL= ', U0*L/NU
		read(*,*)

!----------------- Initial fields -----------------------------

        DO I=1,NI
          DO J=1,NJ
            U_n(I,J)=U0
            V_n(I,J)=1.e-5
            P_n(I,J)=P0
          ENDDO
        ENDDO

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Prandtl equations'      
        call  Prandtl(U_n, V_n, dx, dy, Nu, smax, eps)
		Cf(2:) = tw(mu, U_n(1:,:), dy)*2/(R0*U0*U0)
		Cft(2:) = 0.664/sqrt(U0*R0*X_Node(2:,1)/MU)

 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data' 
        Open(newunit = IO,FILE='data_pr.plt')
        Call Output_Fields(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
        Close(IO)
		
		j = (NI-1)/5
		
		do i = 1,4
		write(*,*)'xplane ', X_Node(j*i+1,1)
		call Blazius(U_n, Y_Node, X_Node, X_Node(j*i+1,1))
		end do
		
		call Output_wall(X_Node(2:,1),Cf(:), Cft(:))
     
        END PROGRAM
   
