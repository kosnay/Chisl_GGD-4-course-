	module InOutput
implicit none
contains
!************************************************************************************************
       SUBROUTINE Output_Fields(IO,NI,NJ,X,Y,U,V,P)
 
         INTEGER NI,NJ,IO, IU
         REAL,DIMENSION(NI,NJ):: X,Y
         REAL,DIMENSION(NI,NJ):: U,V,P

         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI,1:NJ)
         Write(IO,'(100E25.16)') V(1:NI,1:NJ)
         Write(IO,'(100E25.16)') P(1:NI,1:NJ)
		 

       END  SUBROUTINE
	   
	   SUBROUTINE Output_wall(X,Cf,Cft)
 
         INTEGER IO, NI, i
         REAL,intent(in):: X(:), Cf(:), Cft(:)
		 
		 NI = size(CF)
		 
		 open(newunit = IO, file = 'wall.plt')
		 
		 10 format(E16.9E2,2(1x,E16.9E2))

         Write(IO,*) 'VARIABLES = "X", "Cf", "Cft"' 
         Write(IO,*) 'Zone i=',NI,' j=1'
		 do i = 1, NI
			write(IO,10) x(i), Cf(i), Cft(i)
		 end do	
		 close (io)
		 

       END  SUBROUTINE 
end module

!************************************************************************************************
