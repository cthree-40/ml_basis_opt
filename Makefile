# Makefile for nn_basis_opt
#
#
FC := gfortran
SDIR := ./source
BDIR := ./bin

all:
	cd $(SDIR); $(FC) -o compute_rmse_density.x compute_rmse_density.f90
	cd $(SDIR); $(FC) -o compute_rmse_density_3D.x compute_rmse_density_3D.f90
	cd $(SDIR); $(FC) -o cubedata_into_training_data.x cubedata_into_training_data.f90
	cd $(SDIR); $(FC) -o wfdata_into_training_data.x wfdata_into_training_data.f90
	cp $(SDIR)/compute_rmse_density.x $(BDIR)
	cp $(SDIR)/compute_rmse_density_3D.x $(BDIR)
	cp $(SDIR)/cubedata_into_training_data.x $(BDIR)
	cp $(SDIR)/wfdata_into_training_data.x $(BDIR)

clean:
	cd $(SDIR); rm *.o
