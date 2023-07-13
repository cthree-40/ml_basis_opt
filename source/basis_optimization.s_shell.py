#
# create_qchem_file_fhf: Create the qchem input file for FHF
# calculations.
# Input:
#  var = values of parameters
#
def create_qchem_file_fhf(var):
    
    # Find deck file 'fhf.input_deck' in current directory.
    # Append basis set information to file.
    deck_file = open("fhf.input_deck","a")
    deck_file.write("$neo_basis\nH     3\n")
    for i in range(len(var)):
        deck_file.write("S   1   1.0\n")
        deck_file.write(" "+str(var[i])+"   1.0000D+00\n")
    
    deck_file.close()

#
# generate_initial_dataset: Generate an initial set of training
# data. This will be properties from QChem jobs for different
# values of the parameters.
# Input:
#  init_data_size = number of data points
#  num_param      = number of parameters of objective function
#
def generate_initial_dataset(init_data_size, num_param):
    
    # Initialize the values in var_ij. Get set of random values
    # for parameters. Valid values for the parameters occur within
    # bounds.
    var = [[0.0] * num_param] * init_data_size

    for i in range(init_data_size):
        val = objective_function_value(var[i])


#
# objective_function_value: Compute the value of the objective
# function for the parameters.
#
# Info: This objective function computes the RMS error between
# QChem and FGH densities for x and z axes.
#
#  f(x,y,..,z) = SQRT(w_i*(p(x,y,..,z)_i - ref_i)^2)/SQRT(64)
#  w_i         = weight of data at point i
#  p(x,y,..,z) = density at point i
#  ref_i       = FGH density at point i
#
# Input:
#  var = parameters
# Output:
#  val = f(x,y,...,z)
def objective_function_value(var):
    
    # Compute QChem job with values and extract densities
    qchem_job(var)
    
    # Compute RMSE of density using program
    cl('''compute_rmse_density.x fhf_dens.data fhf_dens.fgh.data > fhf.rmse''')
    rmse_file = open("fhf.rmse", "r")
    rmse_file_lines = (rmse_file.read().splitlines())
    val = rmse_file_lines[0]
    return float(val)
    

#
# qchem_job: Compute a QChem job with basis set parameters given by
# input. Return the density. 
# Input:
#  var = basis set parameters
#
def qchem_job(var):
    # Create the input file and run QChem. Then extract the
    # density information, and print to file 'fhf_dens.data'
    create_qchem_file_fhf(var)
    cl('''qchem -nt 24 fhf.input > fhf.output''')
    cl('''cubedata_into_xz_training_data.x den_p_0.cube > fhf_dens.data''')

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
    init_data_size = 100

    # Generate initial data set 
    generate_initial_dataset(init_data_size, num_param)
    
        
    
