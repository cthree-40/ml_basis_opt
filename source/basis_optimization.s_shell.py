import random
import os
import sys
import subprocess
import shlex
import numpy as np


#---------------------------------------------------------------------
# GLOBAL PARAMETERS UNIQUE TO THIS PROGRAM

# This is the directory containing source/ and bin/ and utilities/ for this
# program.
PROGRAM_HOME = "/home/clm96/pi_project/software/basis_opt/malbon_optimizer/nn_basis_opt/"

# Systems used to train model. [Name, # of density pts]
MOLEC_SYS = [["hcn", 64], ["h2o", 96], ["hehhe", 64]]

#---------------------------------------------------------------------





#
# command_line: Issue command to command line
# Input:
#  command = command line as string
# Output:
#  cl_return = returned value from command line
def command_line(command):
    arg = shlex.split(command)
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    print (output)
    return output
    
#
# create_qchem_file: Create a qchem input file
# Input:
#  molec = molecule to create file for
#  var   = parameters
#
def create_qchem_file(molec, var):

    # Create name as string to open input file deck and input file to make.
    # Copy deck file to make new file. Append basis set information to
    # this file.
    deck_fname = molec+".input_deck"
    qcin_fname = molec+".input"
    
    command_line('''cp '''+deck_fname+''' '''+qcin_fname)
    input_file = open(qcin_fname, "a")
    
    input_file.write("$neo_basis\nH    3\n")
    for i in range(len(var)):
        input_file.write("S   1   1.0\n")
        input_file.write(" %.5f   1.0000D+00\n" % var[i])
    input_file.write("****\n$end")

    input_file.close()



#
# generate_initial_dataset: Generate an initial set of training
# data. This will be properties from QChem jobs for different
# values of the parameters.
# Input:
#  init_data_size = number of data points
#  num_param      = number of parameters of objective function
#  run            = execute QChem jobs
#
def generate_initial_dataset(init_data_size, num_param, run):
    
    # We need random numbers.
    random.seed(100001)
    

    # Initialize the values in var_ij. Get set of random values
    # for parameters. Valid values for the parameters occur within
    # bounds.
    v1 = np.random.uniform(25.0, 50.0, size=(init_data_size))
    v2 = np.random.uniform(0.01, 24.9, size=(init_data_size))
    var = np.vstack((v1,v2)).T

    rmse_val = []
    for i in range(init_data_size):

        if (run):
            rmse_val.append(objective_function_value(var[i]))
        else:
            for j in range(len(MOLEC_SYS)):
                create_qchem_file(MOLEC_SYS[j][0], var[i])
                command_line("cp "+MOLEC_SYS[j][0]+".input "+MOLEC_SYS[j][0]+".input."+str(i))
        

    # If all we require is inputs, leave
    if not (run):
        return
        
    # Process and save data
    training_data_file = open("training.dat", "w")
    for i in range(init_data_size):
        training_data_file.write(" %10.5f" % var[i][0])
        training_data_file.write(" %10.5f" % var[i][1])
        training_data_file.write(" %15.8f\n" % rmse_val[i])
    training_data_file.close()

            

#
# objective_function_value: Compute the value of the objective
# function for the parameters.
#
# Info: This objective function calls program to compute RMSE between
# QChem and FGH densities.
#
# Input:
#  var    = parameters
# Output:
#  val = f(x,y,...,z)
def objective_function_value(var):
    
    # Compute QChem job with values and extract densities
    for i in range(len(MOLEC_SYS)):
        qchem_job(var, MOLEC_SYS[i][0])
        
    # Compute RMSE value for each system
    rmse_prog = PROGRAM_HOME+"/bin/compute_rmse_density.x"
    rmse = []
    for i in range(len(MOLEC_SYS)):
        # Name of system and number of points. Stored as local variables for convenience
        name = MOLEC_SYS[i][0]
        npts = MOELC_SYS[i][1]

        command_line(rmse_prog+" "+name+"_dens.data "+name+"_dens.fgh.data "+npts+" > "+name".rmse")
        rmse_file = open(name+".rmse", "r")
        rmse.append(rmse_file.read())
        rmse_file.close()

    # Compute RMSE value
    val = 0.0
    for i in range(len(rmse)):
        val = val + rmse[i]

    return val
    

#
# qchem_job: Compute a QChem job with basis set parameters given by
# input. Return the density. 
# Input:
#  var = basis set parameters
#  molec = moelcule to run Qchem for
def qchem_job(var, molec):
    # Create the input file and run QChem. Then extract the
    # density information, and print to file '*_dens.data'
    create_qchem_file(molec, var)
    command_line('''qchem -nt 24 '''+molec+'''.input > '''+molec+'''.output''')
    cdxz_prog = PROGRAM_HOME+"/bin/cubedata_into_xz_training_data.x"
    command_line(cdxz_prog+''' den_p_0.cube > '''+molec+'''_dens.data''')    

#
# Main Program
#
# 1. Generate initial data set (Y/N)
#    a. Run Qchem
#    b. Compute objective function
#
# 2. 

if __name__ == "__main__":
    

    # Parameters
    num_param = 2
    init_data_size = 20

    # Generate initial data set
    generate_initial_dataset(init_data_size, num_param, False)

    
    
    
