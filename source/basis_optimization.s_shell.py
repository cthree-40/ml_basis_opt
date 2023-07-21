import random
import os
import sys
import time
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

#
# To make this program a bit easier to use, and to not have to multiple copies
# the code is designed to use with a shell scrip that will create the proper
# code using sed.
# From and input file, the variables are printed below.
#
### GLOBAL VAR FROM INPUT ###


# 
# Default parameters commented out below.
#

# Systems used to train model. [Name, # of density pts, cubedata jobtype]
#MOLEC_SYS = [["hcn", 64, 2], ["hehhe", 64, 2]]

# Number of parameters being optimized
#NUM_PARAM = 2

# Orbital type
#ORBITAL_TYPE = "S"

# Build training data set
# True  = Evaluate F(x) for every set of x and build training.dat file
# False = Create input files for every set of x, but do not evaluate
#BUILD_TRAINING_SET = True
#BUILD_TRANING_SET_ONLY = True

# Run GPR fit?
# RUN_GPR = True 

# Initial training data set size
#INIT_DATA_SIZE = 30

# Build testing data set
#BUILD_TESTING_SET = True
#BUILD_TESTING_SET_ONLY = True

# Size of testing data set
#TEST_DATA_SIZE = 3

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
        input_file.write(str(ORBITAL_TYPE)+"   1   1.0\n")
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
    
    # See if training file exists, if it does return.
    if os.path.isfile("training.dat"):
        return

    # We need random numbers.
    random.seed(100001)
    

    # Initialize the values in var_ij. Get set of random values
    # for parameters. Valid values for the parameters occur within
    # bounds.
    var = np.ndarray(shape=(init_data_size, num_param), dtype=float)
    bnds = get_param_bounds(num_param, 50.0, 0.001)
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
    random.seed(200001)

    # Initialize the values in var_ij. Get set of random values
    # for parameters. Valid values for the parameters occur within
    # bounds.
    var = np.ndarray(shape=(test_data_size, num_param), dtype=float)
    bnds = get_param_bounds(num_param, 50.0, 0.001)
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
        testing_data_file.write(" %10.5f" % var[i][0])
        testing_data_file.write(" %10.5f" % var[i][1])
        testing_data_file.write(" %15.8f\n" % rmse_val[i])
    testing_data_file.close()

#
# gp_objfcn: Evaluate the GP objective function
# Input:
#  gp = GP object
#  X  = np.divide(1/X)
def gp_objfcn(X, gp):
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
    var = 1 #I select Radial-basis kernel, user can change it
    i = int(var)
    ks = k_[i]
    print ('You selected --->', i, kt_[i])
    return ks

    
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
# get_param_bounds: Get lower and upper bounds for each parameter.
# Input:
#  num_param = number of parameters
#  hi        = greatest upper bound
#  lo        = least lower bound
# Output:
#  bnds = Upper and lower bounds for each parameter
#
def get_param_bounds(num_param, hi, lo):
    # Compute bound range. Loop over each parameter computing the
    # bounds for each and adding to array.
    #upper_bound = hi
    #bound_range = upper_bound / float(num_param)
    #lower_bound = upper_bound - bound_range
    #bnds = []
    #for i in range(num_param):
    #    bnds.append((lower_bound, upper_bound))
    #    upper_bound = upper_bound - bound_range - 0.01
    #    lower_bound = upper_bound - 2.0 * bound_range
    #    if (i == (num_param - 2)):
    #        lower_bound = lo
    bnds = []
    for i in range(num_param):
        bnds.append((lo, hi))

    return bnds

#
# qchem_job: Compute a QChem job with basis set parameters given by
# input. Return the density. 
# Input:
#  var = basis set parameters
#  molec = moelcule to run Qchem for
#  cjobtype = cubedata_into*x jobetype (1=x, 2=x+z, 3=xyz)
#
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
    ck = C(1.0, (1e-3, 1e3))
    kernels_m = [RBF(length_scale=np.ones(d), length_scale_bounds=(1e-3, 1e+3)),
                 RQ(length_scale=1.0, alpha=0.1,length_scale_bounds=(1e-3, 1e3), alpha_bounds=(1e-3, 1e3)),
                 Matern(length_scale=np.ones(d), length_scale_bounds=(1e-3, 1e+3),nu=2.5)]
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
    gp = GaussianProcessRegressor(kernel = k, n_restarts_optimizer = 12)
    gp.fit(np.atleast_2d(X), Y)
    print("Log-marginal-likelihood (GP): %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))
    
    print("Prediction with GP model")
    y_pred, sigma = gp_prediction(gp, Xt)
    ygp = np.column_stack((Yt, y_pred, sigma))
    dpred = np.column_stack((Xt, ygp))
    print("Save predicted points")
    np.savetxt(gp_test_output, ygp)
    rmse = np.sqrt(mse(Yt, y_pred))
    print("RMSE = ", rmse)
    print(" ----------------------------------\n")
    
    # Global optimization search
    opt_iter = 10
    for j in range(1):
        
        # Get lower/upper bounds for each parameter
        bnds = get_param_bounds(NUM_PARAM, 50.0, 0.001)

        minimizer = {"method": "SLSQP", "args":gp,"bounds":bnds}
        res = basinhopping(gp_objfcn,X[10*(j),:],minimizer_kwargs=minimizer,niter=2000)
        for i in range(NUM_PARAM):
            var[NUM_PARAM*j+i] = res.x[i]
            #X[0,i] = np.divide(1,res.x[i])
            X[0,i] = res.x[i]
            out, sigma = gp_prediction(gp,X[0,:])
            result[j]=out[0]
    for i in range(NUM_PARAM):
        print(" %10.5f" % var[i])
    print(" %15.8f" % result[0])

    print("--- %s seconds ---" % (time.time() - start_time))

        


#
# Main Program
#
# Outline:
#
#
if __name__ == "__main__":
    

    # Generate initial data set
    generate_initial_dataset(INIT_DATA_SIZE, NUM_PARAM, BUILD_TRAINING_SET)
    if (BUILD_TRAINING_SET_ONLY and not BUILD_TESTING_SET_ONLY):
        print("Finished building training set input. Exiting...")
        exit

    # Generate test data set
    generate_testing_dataset(TEST_DATA_SIZE, NUM_PARAM)
    if (BUILD_TESTING_SET_ONLY):
        print("Finished building testing set input. Exiting...")
        exit
    
    if not (RUN_GPR):
        exit

    # Fit GP model and obtain optimized parameters from these results.
    var = [0.0] * NUM_PARAM
    result = [0.0] * 1
    train_gp_and_return_opt(var, result)

    # Check if RMSE is acceptable
    while (abs(result[0]) > 0.0001):
        # Update training and test data sets
        f_training = "training.dat"
        f_testing  = "testing.dat"
        X = np.loadtxt(f_training, usecols=range(NUM_PARAM))
        Y = np.loadtxt(f_training, usecols=(NUM_PARAM,))
        n_points_data, n_parameters_data = X.shape
        for i in range(1):
            # Add new parameters and results
            a = np.array([var])
            X = np.concatenate((X, a), axis=0)
            b = np.array([result[0]])
            Y = np.concatenate((Y, b), axis=0)
        data = np.column_stack((X,Y))
        np.savetxt(f_training, data, fmt=" %10.5f %10.5f %15.8f")
        np.savetxt(f_testing,  data, fmt=" %10.5f %10.5f %15.8f")
        train_gp_and_return_opt(var, result)
    print(var, result[0])
    
    
