import random
import os
import sys
import subprocess
import shlex
import numpy as np

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel,\
	ConstantKernel as C,\
	RationalQuadratic as RQ,\
	ExpSineSquared as ESS,\
	DotProduct as DP
from sklearn.metrics import mean_squared_error as mse 

from scipy.optimize import minimize, rosen, rosen_der,shgo, differential_evolution,basinhopping


#---------------------------------------------------------------------
# GLOBAL PARAMETERS UNIQUE TO THIS PROGRAM

# This is the directory containing source/ and bin/ and utilities/ for this
# program.
PROGRAM_HOME = "/home/clm96/pi_project/software/basis_opt/malbon_optimizer/nn_basis_opt/"

# Systems used to train model. [Name, # of density pts, cubedata jobtype]
MOLEC_SYS = [["hcn", 64, 2], ["hehhe", 64, 2]]

# Number of parameters being optimized
NUM_PARAM = 2

# Initial training data set size
INIT_DATA_SIZE = 30

# Number of OMP Threads
NTHREADS = os.getenv("OMP_NUM_THREADS")
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
# create_nbox_files: Create the nbox_data.txt and nbox_npts.txt files the
# given molecule
# Input:
#  molec = System  
#
def create_nbox_files(molec):
    cmdl_arg = "cp "+molec+".nbox_npts.txt nbox_npts.txt"
    command_line(cmdl_arg)
    cmdl_arg = "cp "+molec+".nbox_data.txt nbox_data.txt"
    command_line(cmdl_arg)


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
    var = np.ndarray(shape=(init_data_size, num_param), dtype=float)
    upper_bound = 50.0
    bound_range = upper_bound / float(num_param)
    lower_bound = upper_bound - bound_range
    print("Generating training set with bounds: ")
    for i in range(num_param):
        print("("+str(lower_bound)+", "+str(upper_bound)+"), ")
        var[:,i] = np.random.uniform(lower_bound, upper_bound, size=(init_data_size))
        # Update lower and upper bounds
        upper_bound = lower_bound - 0.01
        lower_bound = lower_bound - bound_range
        if (i == num_param - 1):
            lower_bound = 0.01 # Make last lower bound
    print("\n")
    sys.stdout.flush()

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
        qchem_job(var, MOLEC_SYS[i][0], MOLEC_SYS[i][2])
        
    # Compute RMSE value for each system
    rmse_prog = PROGRAM_HOME+"/bin/compute_rmse_density.x"
    rmse = []
    for i in range(len(MOLEC_SYS)):
        # Name of system and number of points. Stored as local variables for convenience
        name = MOLEC_SYS[i][0]
        npts = MOLEC_SYS[i][1]

        command_line(rmse_prog+" "+name+"_dens.data "+name+"_dens.fgh.data "+str(npts)+" > "+name+".rmse")
        rmse_file = open(name+".rmse", "r")
        rmse.append(float(rmse_file.read()))
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
#  cjobtype = cubedata_into*x jobetype (1=x, 2=x+z, 3=xyz)
def qchem_job(var, molec, cjobtype):
    # Create the input file and run QChem. Then extract the
    # density information, and print to file '*_dens.data'
    create_qchem_file(molec, var)
    create_nbox_files(molec)
    command_line("qchem -nt "+str(NTHREADS)+" "+molec+".input > "+molec+".output")
    cdxz_prog = PROGRAM_HOME+"/bin/cubedata_into_training_data.x"
    command_line(cdxz_prog+" den_p_0.cube "+str(cjobtype)+" > "+molec+"_dens.data")    
    # Save logs
    command_line("cat "+molec+".output >> qc."+molec+".log")

#
# Main Program
#
# 1. Generate initial data set (Y/N)
#    a. Run Qchem
#    b. Compute objective function
#
# 2. 

if __name__ == "__main__":
    

    # Generate initial data set
    generate_initial_dataset(INIT_DATA_SIZE, NUM_PARAM, True)

    
    
    
