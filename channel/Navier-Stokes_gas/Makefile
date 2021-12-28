FC = ifort
FFLAGS = -O2 -warn all -Wextra -stand f18 -r8 -qopenmp -checkall
SRC = mod_Functions.f90 mod_Output.f90 mod_NavierStokes.f90 Main.f90
OBJ = $(SRC:.f90=.o)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
	
run: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
	
clean:
	@rm -f *.mod *.o run