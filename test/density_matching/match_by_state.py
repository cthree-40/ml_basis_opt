import sys
import math
import numpy as np

if __name__ == "__main__":

    print("Matching states by density.\n")

    dirs = np.loadtxt("hcn.density_matching.txt", usecols=range(1,4))

    nstates = dirs.shape[0]

    for i in range(1,nstates):

        
        if (dirs[i,0] == 0 and dirs[i,1] == 0 and dirs[i,2] == 1):

            print("Comparing Z axis...")
            
            nstart = 64
            nfinal = 95+1
            
            densfgh = np.loadtxt("hcn_dens_"+str(i)+".fgh.dat")

            for j in range(1, 7):

                densabi = np.loadtxt("hcn_dens_"+str(j)+".dat")

                # Compute RMSE
                rmse_vec = densfgh[nstart:nfinal, 1] - densabi[nstart:nfinal, 1]
                rmse_dot = np.dot(rmse_vec, rmse_vec) / float(nfinal - nstart)
                rmse = math.sqrt(rmse_dot)

                print("RMSE with state "+str(j)+" = "+str(rmse))
            
                
            
    print("Done.")
