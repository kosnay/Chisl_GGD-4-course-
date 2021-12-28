FC = gfortran
FFLAGS = -O0 -Wall -Wextra -std=f2018 -fdefault-real-8 -lblas -llapack
SRC = mod_Functions.f90 mod_Output.f90 mod_Prandtl.f90 Main.f90
OBJ = $(SRC:.f90=.o)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
	
run: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
	
clean:
	@rm -f *.mod *.o run