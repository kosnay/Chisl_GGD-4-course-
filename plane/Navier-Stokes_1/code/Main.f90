        Program Pr
		
		use Navier, only : NavierStokes
		use Output, only : OutputFields_Cell, Output_wall
		use Functions, only : tw, Blazius
		
        Implicit none
		
        character(len=100),parameter :: input_file='input_params.nml'
        INTEGER I,J,NI, NJ, NITER, IO
        REAL L,H,U0,MU,NU,R0,P0
        REAL dx,dy,CFL,Uref,EPS
        REAL,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL,ALLOCATABLE :: X_Cell(:,:),Y_Cell(:,:)
        REAL,ALLOCATABLE :: P_c(:,:),U_c(:,:),V_c(:,:)
		REAL,ALLOCATABLE :: Cf(:), Cft(:)
 
        namelist /input_params/ NI,NJ,L,H,U0,MU,R0,P0,EPS,NITER,CFL,Uref
		open(newunit = IO, file = input_file, action = 'read')
		read(IO, nml = input_params)
		close(IO)
		
		write(*,*) 'Read input file - Done'

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! cell centers X-coordinates
        allocate(Y_Cell(0:NI,0:NJ)) ! cell centers Y-coordinates

!------ Cell-centered variables
        allocate(U_c(0:NI,0:NJ))   ! Velocity U
        allocate(V_c(0:NI,0:NJ))   ! Velocity V
        allocate(P_c(0:NI,0:NJ))   ! Pressure
		
		allocate(Cf(0:NI), Cft(0:NI))
		
		Write(*,*) 'Allocation - Done'

        dx=L/(NI-1)
        dy=H/(NJ-1)

!------ Coordinate of nodes
        DO I=1,NI
          DO J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          END DO
        END DO

!------ Coordinate of cell centers
        X_Cell(0,1:NJ)=-dx/2
        Y_Cell(0,1:NJ)=Y_Node(1,1:NJ)+dy/2
        X_Cell(1:NI,0)=X_Node(1:NI,1)+dx/2
        Y_Cell(1:NI,0)=-dy/2
        DO I=1,NI
          DO J=1,NJ
            X_Cell(I,J)=X_Node(I,J)+dx/2
            Y_Cell(I,J)=Y_Node(I,J)+dy/2
          END DO
        END DO

!----------------- Parameters ------------------------

        NU=MU/R0
		
		write(*,*) 'Parameters:'
        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'ReL= ', U0*L/NU
        pause    

!----------------- Initial fields -----------------------------

        DO I=0,NI
          DO J=0,NJ
            U_c(I,J)=U0
            V_c(I,J)=1.0e-5
            P_c(I,J)=0.0
          ENDDO
        ENDDO

!---------------- Solve Navier-Stokes equations ---------------------
 
        write(*,*) 'Solve Navier-Stokes equations' 
        Call NavierStokes(P_c(:,:),U_c(:,:),V_c(:,:), dx, dy, CFL, U0, P0, NU, Uref, NITER, EPS)
                       
 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data cell (Navier-Stokes)' 
        Open(newunit = IO,FILE='data.plt')
        Call OutputFields_Cell(IO,NI,NJ,X_Node,Y_Node,U_c,V_c,P_c)
        Close(IO)
		
		j = (NI-1)/5
		
		do i = 1,4
		write(*,*)'xplane ', X_Node(j*i+1,1)
		call Blazius(U_c, Y_Cell, X_Cell, X_Cell(j*i+1,1))
		end do
		
		Cf = tw (mu, U_c(1:NI-1, 0:1), dy)*2/(R0*(U0*U0))
		Cft = 0.664/sqrt(U0*X_cell(1:NI-1,1)*R0/MU)
		call Output_wall(X_cell(1:NI-1,1),Cf, Cft)
     
        END PROGRAM
   
