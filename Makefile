FC = gfortran -O3

.phony: all clean

all: xproduct.exe invert.exe

xproduct.exe: crossproduct.f90
	@${FC} crossproduct.f90 -o xproduct.exe

invert.exe: matrix-inverse.f90
	@${FC} matrix-inverse.f90 -o invert.exe

clean: 
	@rm *.exe