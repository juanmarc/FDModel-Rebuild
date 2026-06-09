#build script for FD model
rm *.o
gfortran -c -fcheck=all basic.state.f
gfortran -c -fcheck=all dncarpois.f
gfortran -c -fcheck=all grid.initialize.f
gfortran -c -fcheck=all modelOutput.f
gfortran -c -fcheck=all pois.solver.f
gfortran -c -fcheck=all poisson.setup.f
gfortran -c -fcheck=all setup.f
gfortran -fcheck=all fd.driver.f -o fdtest basic.state.o dncarpois.o \
	grid.initialize.o modelOutput.o pois.solver.o poisson.setup.o setup.o

