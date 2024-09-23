import math
import random
import os
import sys
import time
import subprocess
import shlex
import numpy as np

from sklearn import preprocessing
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel,\
	ConstantKernel as C,\
	RationalQuadratic as RQ,\
	ExpSineSquared as ESS,\
	DotProduct as DP
from sklearn.metrics import mean_squared_error as mse 

from scipy.optimize import minimize, rosen, rosen_der,shgo, differential_evolution,basinhopping


#---------------------------------------------------------------------


#
# To make this program a bit easier to use, and to not have to multiple copies
# the code is designed to use with a shell scrip that will create the proper
# code using sed.
# From and input file, the variables are printed below.
#


#
# DEFAULT PARAMETERS ##############################################
#

# This is the directory containing source/ and bin/ and utilities/ for this
# program.
PROGRAM_HOME = "/home/clm96/pi_project/software/basis_opt/malbon_optimizer/ml_basis_opt"

# Jobtype for program execution
#  gen_training: Generate training input files for batch submission
#  gen_testing:  Generate testing input files for batch submission
#  run_gpr:      Run GPR
#  all:          Do all three
JOBTYPE = "run_gpr"

# Build the surface, or just return parameters
BUILD_SURFACE = True

# Maximum number of iterations
MAX_ITER = 10

# Random seeds
TRAIN_RSEED = 100001
TEST_RSEED  = 200001

# Kernel to use
# 0 = Matern
# 1 = RBF
# 2 = RQ
MLKERNEL = 1

# Systems used to train model. 
#  [Name, # of density pts, cubedata jobtype, # states]
MOLEC_SYS = [["hcn", 64, 2, 1], ["hehhe", 64, 2, 1]]

# Number of parameters being optimized
NUM_PARAM = 2

# Bounds for parameters
PBOUNDS = [(0.01, 20.0), (0.01, 50.0)]

# Orbital type
ORBITAL_TYPE = ["S","S"]

# Penalty function weights
PENFCN_WEIGHTS = {"density" : 10.0, "excited states" : 1.0, "ground state" : 1.0, "states" : 0.0}

# Build training data set
# True  = Evaluate F(x) for every set of x and build training.dat file
# False = Create input files for every set of x, but do not evaluate
BUILD_TRAINING_SET = True

# Initial training data set size
INIT_DATA_SIZE = 30

# Build testing data set
BUILD_TESTING_SET = True

# Size of testing data set
TEST_DATA_SIZE = 3

# Compute penalty function in parallel
PARALLEL = "n"

# Use inverse of X (1/X)
USE_XINV = False

# Enforce linear independence by ensuring coefficients are separated by
# linind_coeftol amount. Sets P(x) = LININD_PFCNVAL 
ENFORCE_LININD = False
LININD_COEFTOL = 0.5
LININD_PFCNVAL = 30000

# Do not accept as predicted minima P(x) < 0.0 (i.e. holes) that may be
# at edge of range.
NOEDGE_MINIMA = False

# Check P(x) function
CHECK_PENFCN = False

# Penalty function tolerance for convergence
PF_TOL = 1.0

# Global optimization searches
NMIN_SEARCH = 5

# Start one search from minimum of data set
START_FROM_MIN = False

# Scaling of P(x)
YSCALE = 1.0

# String for beginning of basis
BASBEGINSTR="$neo_basis\nH    3\n"

# Build even tempered basis
EVEN_TEMPERED_BASIS = False
ET_BASIS_NUM_FCNS = [5, 4, 3, 2, 1]

#######################################################################
### GLOBAL VAR FROM INPUT ###


# Number of OMP Threads
NTHREADS = os.getenv("OMP_NUM_THREADS")
#---------------------------------------------------------------------


# check_linind: Check linear independence of functions
# Returns:
#  False = functions could be linearly independent
#  True  = functions could be linearly dependent
def check_linind(X):

    if ENFORCE_LININD == False:
        return False
    
    # Loop over paramters. If parameters are of same species check
    # if paramters are separated by at least LININD_COEFTOL.
    for i in range(NUM_PARAM - 1):
        if ORBITAL_TYPE[i] == ORBITAL_TYPE[i+1]:
            if abs(X[i] - X[i+1]) < LININD_COEFTOL:
                return True
    
    return False

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
    
    input_file.write(BASBEGINSTR)

    if EVEN_TEMPERED_BASIS == False:

        for i in range(len(var)):
            input_file.write(str(ORBITAL_TYPE[i])+"   1   1.0\n")
            input_file.write(" %.5f   1.0000D+00\n" % var[i])

    else:

        vix = 0 # start var indexing, we will jump by two for each a,b pair
        for i in range(len(ET_BASIS_NUM_FCNS)):
            for j in range(ET_BASIS_NUM_FCNS[i]):
                expcoef = var[vix]*math.pow(var[vix+1],j) # j = n-1 = 0 .. n-1
                input_file.write(str(ORBITAL_TYPE[i])+"   1   1.0\n")
                input_file.write(" %.5f   1.0000D+00\n" % expcoef)

            vix = vix + 2 # increment to next a,b pair


    input_file.write("****\n$end\n")        
    input_file.close()

#
# create_pjob_shell_script: Create a shell script to execute jobs in
# parallel
#
def create_pjob_shell_script(fname):


    shell_file = open(fname, "w")

    PCTHREADS = int(NTHREADS) // 2 # Only two jobs in parallel
    
    shell_file.write("#!/bin/bash\n")

    # First system
    mname = MOLEC_SYS[0][0]
    shell_file.write("cd "+mname+"\n")
    shell_file.write("qchem -nt "+str(NTHREADS)+" "+mname+".input > "+mname+".output\n")
    shell_file.write("cd ../\n")

    # Parallel systems
    mname = MOLEC_SYS[1][0]
    shell_file.write("cd "+mname+"\n")
    shell_file.write("qchem -nt "+str(PCTHREADS)+" "+mname+".input > "+mname+".output &\n")
    shell_file.write("cd ../\n")
    mname = MOLEC_SYS[2][0]
    shell_file.write("cd "+mname+"\n")
    shell_file.write("qchem -nt "+str(PCTHREADS)+" "+mname+".input > "+mname+".output &\n")
    shell_file.write("cd ../\n")

    shell_file.write("wait\n")

    # After jobs finish, extract density information and copy output files to parent directory
    for i in range(len(MOLEC_SYS)):
        mname = MOLEC_SYS[i][0]
        cjtype= MOLEC_SYS[i][2]
        shell_file.write("cd "+mname+"\n")

        # Compute density and copy to parent directory
        cdxz_prog = PROGRAM_HOME+"/bin/cubedata_into_training_data.x"
        cdxz_cmdl = cdxz_prog+" den_p_0.cube "+str(cjtype)+" > "+mname+"_dens.data"    
        shell_file.write(cdxz_cmdl+"\n")
        shell_file.write("cp "+mname+"_dens.data ../\n")
        
        # Copy output file to parent directory
        shell_file.write("cp "+mname+".output ../\n")

        # Leave directory, returning to parent directory
        shell_file.write("cd ../\n")
        
    shell_file.close()

    # Make executable
    command_line("chmod +x "+fname)
    
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
    
    # See if training file exists, if it does return.
    if os.path.isfile("training.dat"):
        print("Training data file found! Aborting...\n")
        return

    # We need random numbers.
    np.random.seed(TRAIN_RSEED)

    # Initialize the values in var_ij. Get set of random values
    # for parameters. Valid values for the parameters occur within
    # bounds.
    var = np.ndarray(shape=(init_data_size, num_param), dtype=float)
    bnds = get_param_bounds()
    print("Generating training set with bounds: ")
    for i in range(num_param):
        hi = bnds[i][1]
        lo = bnds[i][0]
        print("("+str(lo)+", "+str(hi)+"), ")
        var[:,i] = np.random.uniform(lo, hi, size=(init_data_size))

    print("\n")
    sys.stdout.flush()

    rmse_val = []
    for i in range(init_data_size):

        if (run):
            rmse_val.append(objective_function_value(var[i]))
        else:
            # save var[i] info
            save_var_to_indexed_file(var[i], i)
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
# generate_testing_dataset: Generate test set of points to evaluate model
# Input:
#  test_data_size = number of test points
#  num_param      = number of parameters of objective function
#
def generate_testing_dataset(test_data_size, num_param):
    
    # See if testing data exists. If it does, return
    if os.path.isfile("testing.dat"):
        return

    # We need random numbers. (Different seed from training data set)
    np.random.seed(TEST_RSEED)

    # Initialize the values in var_ij. Get set of random values
    # for parameters. Valid values for the parameters occur within
    # bounds.
    var = np.ndarray(shape=(test_data_size, num_param), dtype=float)
    bnds = get_param_bounds()
    print("Generating testing data set with bounds: ")
    for i in range(num_param):
        hi = bnds[i][1]
        lo = bnds[i][0]
        print("("+str(lo)+", "+str(hi)+"), ")
        var[:,i] = np.random.uniform(lo, hi, size=(test_data_size))

    print("\n")
    sys.stdout.flush()
    
    # Evaluate function at these points
    rmse_val = []
    for i in range(test_data_size):
            rmse_val.append(objective_function_value(var[i]))

    # Process and save data
    testing_data_file = open("testing.dat", "w")
    for i in range(test_data_size):
        for j in range(num_param):
            testing_data_file.write(" %10.5f" % var[i][j])
            
        testing_data_file.write(" %15.8f\n" % rmse_val[i])
    testing_data_file.close()

#
# gp_objfcn: Evaluate the GP objective function
# Input:
#  gp = GP object
#  X  = parameters
def gp_objfcn(X, gp):
    if (check_linind(X)):
        y_pred = LININD_PFCNVAL
    if (USE_XINV):
        Xinv = np.divide(1, X)
    y_pred, sigma = gp.predict(np.atleast_2d(X), return_std=True)
    return y_pred

#
# gp_prediction: Return prediction from GP
# Input:
#  gp = trained GP model
#  X  = points to predict
def gp_prediction(gp, X):
    y_pred, sigma = gp.predict(np.atleast_2d(X), return_std=True)
    return y_pred, sigma

#
# kernel_user: User selects kernel
#
def kernel_user():
    k_ = ['mat','rbf','rq']
    kt_ = ['Matern kernel', 'Radial-basis function kernel', 'Rational Quadratic kernel']
    print ('Select kernel type')
    print ('Select kernel type')
    print ('For Mattern function, enter 0')
    print ('For Radial-basis function kernel, enter 1')
    print ('For Rational Quadratic kernel, enter 2')
    var = MLKERNEL
    i = int(var)
    ks = k_[i]
    print ('You selected --->', i, kt_[i])
    return ks

    
#
# objective_function_value: Compute the objective function for the parameters
#
def objective_function_value(var):

    # Compute QChem job with values and extract densities and states
    if (PARALLEL == "n"):
        for i in range(len(MOLEC_SYS)):
            qchem_job(var, MOLEC_SYS[i][0], MOLEC_SYS[i][2], MOLEC_SYS[i][3])
    else:
        qchem_pjobs(var)
    
    wt = PENFCN_WEIGHTS
    
    val = 0.0
    if (wt["density"] > 0.0):
        val = val + wt["density"] * objective_function_value_density(var)

    val = val + wt["excited states"] * objective_function_value_excited_states(var)
    val = val + wt["ground state"] * objective_function_value_zeropoint(var)
    val = val + wt["states"] * objective_function_value_all_states(var)
    
    return val

#
# objective_function_value_all_states: Compute the all states component of the
# objective function.
#
# Info: This objective function calls programs to compute RMSE between QChem and 
# FGH energies.
#
# Input:
#  var = parameters
# Output:
#  RMSE for all states
#
def objective_function_value_all_states(var):

    if (PENFCN_WEIGHTS["states"] == 0.0):
        val = 0.0
        return val

    
    # Compute RMSE value for each system
    rmse = []
    for i in range(len(MOLEC_SYS)):
        # Name of system and number of states.
        name = MOLEC_SYS[i][0]
        nsts = MOLEC_SYS[i][3]
        
        # Read in excited states. Arrive converted in au
        states = get_all_states_from_file(name+"_states.data", nsts)
        # Read in reference excited states. Arrive converted to au
        ref_st = get_all_states_from_file(name+"_states.fgh.data", nsts)

        # Compute RMSE
        val = 0.0
        for i in range(nsts):
            
            diff = (states[i] - ref_st[i]) * 219474.63 # Convert to cm-1
            
            val = val + diff * diff

        val = val / float(nsts)
        val = math.sqrt(val)
        rmse.append(val)
        
    # Sum errors
    rmse_val = 0.0
    for i in range(len(MOLEC_SYS)):
        rmse_val = rmse_val + rmse[i]

    return rmse_val


#
# objective_function_value_density: Compute the density component of the objective
# function for the parameters.
#
# Info: This objective function calls program to compute RMSE between
# QChem and FGH densities.
#
# Input:
#  var    = parameters
# Output:
#  val = f(x,y,...,z)
def objective_function_value_density(var):

    if (PENFCN_WEIGHTS["density"] == 0.0):
        val = 0.0
        return val
    
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
# objective_function_value_excited_states: Compute the excited stats component of the
# objective funciton.
#
# Info: This objective function calls programs to compute RMSE between QChem and
# FGH energies.
#
# Input:
#  var  = parameters
# output:
#  RMSE for excited states
#
def objective_function_value_excited_states(var):

    if (PENFCN_WEIGHTS["excited states"] == 0.0):
        val = 0.0
        return val

    
    # Compute RMSE value for each system
    rmse = []
    for i in range(len(MOLEC_SYS)):
        # Name of system and number of states.
        name = MOLEC_SYS[i][0]
        nsts = MOLEC_SYS[i][3]
        nxst = nsts - 1

        if (nxst == 0):
            continue
        
        # Read in excited states. Arrive converted to cm-1
        xstates = get_excited_states_from_file(name+"_states.data", nxst)
        # Read in reference excited states. Arrive converted to cm-1
        ref_xst = get_excited_states_from_file(name+"_states.fgh.data", nxst)

        # Compute RMSE
        val = 0.0
        for i in range(nxst):
            
            diff = xstates[i] - ref_xst[i]
            
            val = val + diff * diff

        val = val / float(nxst)
        val = math.sqrt(val)
        rmse.append(val)
        
    # Sum errors
    rmse_val = 0.0
    for i in range(len(rmse)):
        rmse_val = rmse_val + rmse[i]

    return rmse_val

#
# objective_function_value_zeropoint: Compute the RMSE for zero point 
# energy.
# Input:
#  var = parameter values
# Output:
#  RMSE for zero point energy
#
def objective_function_value_zeropoint(var):

    if (PENFCN_WEIGHTS["ground state"] == 0.0):
        val = 0.0
        return val

    
    rmse = []
    # Compute RMSE value fo each system
    for i in range(len(MOLEC_SYS)):
        
        # Name of system
        name = MOLEC_SYS[i][0]
        
        # Read in ground states. Arrive in au
        e0_qch = get_ground_state_from_file(name+"_states.data")
        # Read in reference ground states. Arrive in au
        e0_ref = get_ground_state_from_file(name+"_states.fgh.data")

        # Get difference in au, and compute RMSE
        val = (e0_qch - e0_ref)
        val = val * val
        val = math.sqrt(val)
        rmse.append(val)
        
    # Sum errors
    rmse_val = 0.0
    for i in range(len(MOLEC_SYS)):
        rmse_val = rmse_val + rmse[i]

    return rmse_val

#
# get_all_states_from_file: Get state energies from
# a file.
# Input:
#  flname = File name containing state information
#  nsts   = Number of excited states
# Output:
#  states = Excited states
#
def get_all_states_from_file(flname, nsts):
    
    states = []

    if os.path.isfile(flname):
        
        efile = open(flname, "r")
        ef_lines = (efile.read().splitlines())
        efile.close()

        nlines = len(ef_lines)
        # Ensure we have all states required
        if (nlines < nsts):
            print("Missing state information.")
            sys.exit("Error message")
            

        for i in range(0, nlines):
            states.append(float(ef_lines[i]))
        
    else:
        
        print("Could not find energy file! "+flname)
        sys.exit("Error message")
    
    return states


#
# get_excited_states_from_file: Get excited state energies from
# a file. Convert these energies to cm-1
# Input:
#  flname = File name containing state information
#  nsts   = Number of excited states
# Output:
#  xstates = Excited states
#
def get_excited_states_from_file(flname, nsts):
    
    xstates = []

    if os.path.isfile(flname):
        
        efile = open(flname, "r")
        ef_lines = (efile.read().splitlines())
        efile.close()

        nlines = len(ef_lines)
        # Ensure we have all states required
        if (nlines <= nsts):
            print("Missing excited state information.")
            sys.exit("Error message")
            

        for i in range(1, nlines):
            xstates.append(float(ef_lines[i]))
        
        # Convert to eV
        for i in range(len(xstates)):
            xstates[i] = (xstates[i] - float(ef_lines[0]))*219474.63
            
    else:
        
        print("Could not find energy file! "+flname)
        sys.exit("Error message")
    
    return xstates

#
# get_ground_state_from_file: Get ground state energies from a file.
# Energy returned is absolute in hartree.
# Input:
#  flname = File name containing state information
# Output:
#  e0 = Ground state energy
#
def get_ground_state_from_file(flname):
    
    if os.path.isfile(flname):
        
        efile = open(flname, "r")
        ef_lines = (efile.read().splitlines())
        efile.close()

        e0 = float(ef_lines[0])

    else:
        
        print("Could not find energy file! "+flname)
        sys.exit("Error message")
    
    return e0

#
# get_param_bounds: Get lower and upper bounds for each parameter.
#
def get_param_bounds():
    return PBOUNDS

#
# qchem_job: Compute a QChem job with basis set parameters given by
# input. 
#
# This function also handles processing file for density / excited states
#
# Input:
#  var = basis set parameters
#  molec = moelcule to run Qchem for
#  cjobtype = cubedata_into*x jobetype (1=x, 2=x+z, 3=xyz)
#
def qchem_job(var, molec, cjobtype, nstates):
    # Create the input file and run QChem. Then extract the
    # density information, and print to file '*_dens.data'
    # After this, get energy information and save it to '*_states.data'
    create_qchem_file(molec, var)
    create_nbox_files(molec)

    command_line("qchem -nt "+str(NTHREADS)+" "+molec+".input > "+molec+".output")
    cdxz_prog = PROGRAM_HOME+"/bin/cubedata_into_training_data.x"
    command_line(cdxz_prog+" den_p_0.cube "+str(cjobtype)+" > "+molec+"_dens.data")    


    qfile = open(molec+".output", "r")
    qf_lines = (qfile.read().splitlines())
    qfile.close()
    nlines = len(qf_lines)
    estart_line = 0
    for i in range(nlines):
        # Loop through qchem output file to find the Final CI Energy
        if "==== Final CI Energy ====" in str(qf_lines[i]):
            estart_line = i
            break
    estart_line = estart_line + 1 # Adjust for first entry
    # Now read data for each state
    energy = []
    for i in range(nstates):
        e_line = qf_lines[estart_line + i]
        e_line = e_line.replace(" CI Energy (au) Root #   "+str(i),"").strip()
        energy.append(e_line)

    efile = open(molec+"_states.data", "w")
    for i in range(nstates):
        efile.write("%s\n" % energy[i])
    efile.close()

    # Save logs
    command_line("cat "+molec+".output >> qc."+molec+".log")

#
# qchem_pjobs: Run QChem jobs in parallel (for last two systems)
#
# Input:
#  var = basis set parameters
#
def qchem_pjobs(var):

    # Create shell script
    fname = "pjob_submit.sh"
    create_pjob_shell_script(fname)

    # Create directories for each molecular system
    for i in range(len(MOLEC_SYS)):
        create_qchem_file(MOLEC_SYS[i][0], var)
        create_nbox_files(MOLEC_SYS[i][0])

        # Copy files into directory
        command_line("mkdir -p "+MOLEC_SYS[i][0])
        command_line("cp "+MOLEC_SYS[i][0]+".input "+MOLEC_SYS[i][0]+"/")
        command_line("cp nbox_*.txt "+MOLEC_SYS[i][0]+"/")

    # Execute pjob shell script
    command_line("./pjob_submit.sh")

    # Get data from each molecule
    for i in range(len(MOLEC_SYS)):
        molec = MOLEC_SYS[i][0]
        qfile = open(molec+".output", "r")
        qf_lines = (qfile.read().splitlines())
        qfile.close()
        nlines = len(qf_lines)
        estart_line = 0
        for j in range(nlines):
            # Loop through qchem output file to find the Final CI Energy
            if "==== Final CI Energy ====" in str(qf_lines[j]):
                estart_line = j
                break

        estart_line = estart_line + 1 # Adjust for first entry
        # Now read data for each state
        energy = []
        nstates = MOLEC_SYS[i][3]
        for j in range(nstates):
            e_line = qf_lines[estart_line + j]
            e_line = e_line.replace(" CI Energy (au) Root #   "+str(j),"").strip()
            energy.append(e_line)
                
        efile = open(molec+"_states.data", "w")
        for j in range(nstates):
            efile.write("%s\n" % energy[j])

        efile.close()

        # Save logs
        command_line("cat "+molec+".output >> qc."+molec+".log")
        

#
# save_var_to_indexed_file: Save the current coefficients to a file.
# This routine is used for processing inputs that were submitted as batch jobs
# and not computed during program execution.
# Input:
#  var = parameters
#    i = file index
#
def save_var_to_indexed_file(var, i):
    f = open("var."+str(i), "w")
    for j in range(NUM_PARAM):
        f.write(" %.5f\n" % var[j])
    f.close()

#
# train_gp_and_return_opt: Train the GP model and return optimized parameters and
# GP results.
# Output:
#  var = parameters
#  result = result of GP
#
def train_gp_and_return_opt(var, result):
    
    start_time = time.time()
    
    training_data_file = "training.dat"
    testing_data_file  = "testing.dat"
    gp_test_output     = "test.output"
    
    # Ensure training and testing files are present
    if not os.path.isfile(training_data_file):
        print("ERROR: Missing training data file!\n")
        sys.exit("Error message")
    if not os.path.isfile(testing_data_file):
        print("ERROR: Missing testing data file!\n")
        sys.exit("Error message")
        
    # Read in training points and testing points
    X  = np.loadtxt(training_data_file, usecols=range(NUM_PARAM))
    Xt = np.loadtxt(testing_data_file,  usecols=range(NUM_PARAM))
    
    Y  = np.loadtxt(training_data_file, usecols=(NUM_PARAM, ))
    Yt = np.loadtxt(testing_data_file,  usecols=(NUM_PARAM, ))
    # Scale
    Y = np.multiply(Y, YSCALE)
    Yt = np.multiply(Yt, YSCALE)
    
    Xinv  = np.divide(1, X)
    Xtinv = np.divide(1, Xt)
    
    XMor  = np.exp(-X  /10.0)
    XtMor = np.exp(-Xt /10.0) 
    
    npoints_training = X.shape[0]
    npoints_testing  = Xt.shape[0]

    print("--------------------------------------------\n")
    print("Training data:\n")
    print(" Number of parameters (dimensions): %10d\n" % NUM_PARAM)
    print(" Number of training points: %10d\n" % npoints_training)
    print("Testing data:\n")
    print(" Number of testing points: %10d\n" % npoints_testing)
    print("--------------------------------------------\n")

    print("\n Total time: %s seconds\n" % (time.time() - start_time))

    # Kernel selection
    d = NUM_PARAM
    ck = C(1.0, (1e-5, 1e+5))
    kernels_m = [RBF(length_scale=np.ones(d), length_scale_bounds=(1e-5, 1e+5)),
                 RQ(length_scale=1.0, alpha=0.1,length_scale_bounds=(1e-5, 1e+5), alpha_bounds=(1e-5, 1e+5)),
                 Matern(length_scale=np.ones(d), length_scale_bounds=(1e-5, 1e+5),nu=2.5)]
    ks = kernel_user()
    if ks == "mat":
        k = ck * kernels_m[2]
    elif ks == "rbf":
        k = ck * kernels_m[0]
    elif ks == "rq":
        k = ck * kernels_m[1]
    wk_bool = "no"
    if wk_bool == "yes" or wk_bool == "YES":
        k = k + ck * WhiteKernel(0.1)

    # Train model
    optdat_file = open("opt.dat", "w")
    result_file = open("result",  "w")
    print("Training GP model")
    gp = GaussianProcessRegressor(kernel = k, n_restarts_optimizer = 12, normalize_y = True)
    if not (USE_XINV):
        gp.fit(np.atleast_2d(X), Y)
    else:
        gp.fit(np.atleast_2d(Xinv), Y)
    
    print("Log-marginal-likelihood (GP): %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

    f_coef = "coeff.dat"
    np.savetxt(f_coef,gp.kernel_.theta,fmt="%.16f")
    np.savetxt("k.dat",gp.L_,fmt="%.16f")
    np.savetxt("alpha.dat",gp.alpha_,fmt="%.16f")

    
    print("Prediction with GP model")
    if not (USE_XINV):
        y_pred, sigma = gp_prediction(gp, Xt)
    else:
        y_pred, sigma = gp_prediction(gp, Xtinv)
    
    ygp = np.column_stack((Yt, y_pred, sigma))
    dpred = np.column_stack((Xt, ygp))
    print("Save predicted points")
    np.savetxt(gp_test_output, ygp)
    rmse = np.sqrt(mse(Yt, y_pred))
    print("RMSE = ", rmse)
    print(" ----------------------------------\n")
    sys.stdout.flush()

    # Optimization
    nsearch = NMIN_SEARCH # Number of guess searches
    var_loc = [0.0] * NUM_PARAM * nsearch
    res_loc = [0.0] * nsearch
    xrows, xcols = X.shape
    # Starting points
    sg = np.random.randint(0, xrows, nsearch)
    datamin = np.argmin(Y)
    print(" Minimum in training data: "+str(datamin)+" , "+str(Y[datamin]))
    if (START_FROM_MIN):
        sg[0] = datamin
        
    for j in range(nsearch):
        
        # Get lower/upper bounds for each parameter
        bnds = get_param_bounds()

        minimizer = {"method": "SLSQP", "args":gp,"bounds":bnds}
        res = basinhopping(gp_objfcn,X[sg[j],:],minimizer_kwargs=minimizer,niter=1000)
        for i in range(NUM_PARAM):
            var_loc[NUM_PARAM*j+i] = res.x[i]
            if (USE_XINV):
                X[0,i] = np.divide(1,res.x[i])
            else:
                X[0,i] = res.x[i]
                
            out, sigma = gp_prediction(gp,X[0,:])
            res_loc[j]=out[0]

    print(" Starting guesses:")
    for j in range(nsearch):
        print(" %10.8f" % Y[sg[j]], end="")
    print("\n")
    print(" Minima:")
    for j in range(nsearch):
        print(" %10.8f" % res_loc[j], end="")
    print("\n")

    # Find lowest minimum
    min_val = 1000000.0
    min_idx = 0
    for j in range(nsearch):
        if (res_loc[j] < min_val):
            # if we are rejecting P(x) = 0.0 edge cases
            if (NOEDGE_MINIMA and res_loc[j] <= 0.0001):
                continue
            
            min_idx = j
            min_val = res_loc[j]
    
    # Print this minimum
    for i in range(NUM_PARAM):
        print(" %10.5f" % var_loc[(NUM_PARAM * min_idx) + i])
    print(" %15.8f" % res_loc[min_idx])
    
    # Return this minimum
    var[0:NUM_PARAM] = var_loc[(NUM_PARAM * min_idx):(NUM_PARAM * min_idx) + NUM_PARAM]
    result[0] = res_loc[min_idx] * YSCALE
        
    print("--- %s seconds ---" % (time.time() - start_time))

        


#
# Main Program
#
# Outline:
#
#
if __name__ == "__main__":
    

    # Generate initial data set
    if (JOBTYPE == "gen_training" or JOBTYPE == "all"):
        generate_initial_dataset(INIT_DATA_SIZE, NUM_PARAM, BUILD_TRAINING_SET)
        
        
    # Generate test data set
    if (JOBTYPE == "gen_testing" or JOBTYPE == "all" or JOBTYPE == "run_gpr"):
        generate_testing_dataset(TEST_DATA_SIZE, NUM_PARAM)
    

    # Fit GP model and obtain optimized parameters from these results.
    if (JOBTYPE == "run_gpr" or JOBTYPE == "all"):
        var = [0.0] * NUM_PARAM
        result = [0.0] * 1
        rmse = [0.0] * 1

        # Check penalty function evaluation
        if (CHECK_PENFCN):
            X_test = np.loadtxt("training.dat", usecols=range(NUM_PARAM))
            Y_test = np.loadtxt("training.dat", usecols=(NUM_PARAM,))
            pt = np.random.randint(0, Y_test.shape)
            test_val = [0.0] * 1
            test_var = X_test[pt,0:NUM_PARAM].tolist()
            test_val[0] = objective_function_value(test_var[0])
            print("Penalty function test at point "+str(pt))
            print("Parameters: "+str(test_var[0]))
            print("  P(x) = %15.8f\n" % test_val[0])
            print("  T(x) = %15.8f\n" % Y_test[pt])
            print("\n")
         
        
        gpr_reliable = False
        iter = 1
        while (iter < MAX_ITER):
            
            # Train GP and return minimzation results and RMSE.
            # If the RMSE is good, this is our minimum. Leave loop.
            train_gp_and_return_opt(var, result)

            # Save the GP predicted minimum, and compute the objective
            # function at this point. Compute difference and print to
            # stdout.
            gpmin = result[0]
            if (BUILD_SURFACE):
                result[0] = objective_function_value(var)
                diff = gpmin - result[0]
                print("GP predicted minimum = %15.8f\n" % gpmin)
                print("Computed minimum     = %15.8f\n" % result[0])
                print("Difference           = %15.8f\n" % diff)
                if ((abs(diff) < PF_TOL) and (gpmin < 100.0)):
                    gpr_reliable = True
                    break
            else:
                # Create new input files
                for i in range(len(MOLEC_SYS)):
                    create_qchem_file(MOLEC_SYS[i][0], var)
                
                break
            
            # Update training and testing data sets. (They are the same)
            f_training = "training.dat"
            f_testing  = "testing.dat"
            X = np.loadtxt(f_training, usecols=range(NUM_PARAM))
            Y = np.loadtxt(f_training, usecols=(NUM_PARAM,))
            n_points_data, n_parameters_data = X.shape
            # Add new parameters and results
            a = np.array([var])
            X = np.concatenate((X, a), axis=0)
            b = np.array([result[0]])
            Y = np.concatenate((Y, b), axis=0)
            data = np.column_stack((X,Y))
            # Build format string
            fmtstr = ""
            for v in range(NUM_PARAM):
                fmtstr = fmtstr+" %10.5f"
            
            fmtstr = fmtstr+" %10.8f"
            # save old training and testing data files
            command_line("cp training.dat training.dat_prev")
            #command_line("cp testing.dat testing.dat_prev")
            # save results
            np.savetxt(f_training, data, fmt=fmtstr)
            #np.savetxt(f_testing,  data, fmt=fmtstr)
            
            # Increment
            iter = iter + 1
        
        if (gpr_reliable):
            print("Valid minimum found!")
            print("Basis:")
            for i in range(NUM_PARAM):
                print(" %10.5f" % var[i], end="")
            print(" %10.8f\n" % result[0])
        else:
            print("Need more points.\n")

    print("Job complete.")
