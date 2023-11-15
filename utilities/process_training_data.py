import sys, getopt
import numpy

def check_for_lindep(data, ns, np, nd, nf, ng, nh):

    nrows = data.shape[0]
    ncols = data.shape[1]

    maxval = data[numpy.argmax(data[:,ncols-1]),ncols-1]

    print("Maximum y-value in data = "+str(maxval))
    
    for pt in range(nrows):

        idx=0
        for i in range(ns - 1):
            if abs(data[pt,idx+i] - data[pt,idx+i+1]) < 0.5:
                data[pt,ncols-1] = maxval
        idx=ns
        for i in range(np - 1):
            if abs(data[pt,idx+i] - data[pt,idx+i+1]) < 0.5:
                data[pt,ncols-1] = maxval
        idx=ns+np
        for i in range(nd - 1):
            if abs(data[pt,idx+i] - data[pt,idx+i+1]) < 0.5:
                data[pt,ncols-1] = maxval
        idx=ns+np+nd
        for i in range(nf - 1):
            if abs(data[pt,idx+i] - data[pt,idx+i+1]) < 0.5:
                data[pt,ncols-1] = maxval
        idx=ns+np+nd+nf
        for i in range(ng - 1):
            if abs(data[pt,idx+i] - data[pt,idx+i+1]) < 0.5:
                data[pt,ncols-1] = maxval
        idx=ns+np+nd+nf+ng
        for i in range(nh - 1):
            if abs(data[pt,idx+i] - data[pt,idx+i+1]) < 0.5:
                data[pt,ncols-1] = maxval
        idx=ns+np+nd+nf+ng+nh
        
def get_arguments_from_cmdl(argv):
    sfcns = 0
    pfcns = 0
    dfcns = 0
    ffcns = 0
    gfcns = 0
    hfcns = 0
    opts, args = getopt.getopt(argv,"s:p:d:f:g:h:")
    for opt, arg in opts:
        if (opt == "-s"):
            sfcns = arg
        elif (opt == "-p"):
            pfcns = arg
        elif (opt == "-d"):
            dfcns = arg
        elif (opt == "-f"):
            ffcns = arg
        elif (opt == "-g"):
            gfcns = arg
        elif (opt == "-h"):
            hfcns = arg
        else:
            print ("Unknown option")
            sys.exit()

    return (int(sfcns), int(pfcns), int(dfcns),
            int(ffcns), int(gfcns), int(hfcns))

def sort_parameters(x, ns, np, nd, nf, ng, nh):
    if (ns > 1):
        start=0
        final=ns
        x[start:final] = numpy.sort(x[start:final])
    if (np > 1):
        start=ns
        final=ns+np
        x[start:final] = numpy.sort(x[start:final])
    if (nd > 1):
        start=ns+np
        final=ns+np+nd
        x[start:final] = numpy.sort(x[start:final])
    if (nf > 1):
        start=ns+np+nd
        final=ns+np+nd+nf
        x[start:final] = numpy.sort(x[start:final])
    if (ng > 1):
        start=ns+np+nd+nf
        final=ns+np+nd+nf+ng
        x[start:final] = numpy.sort(x[start:final])
    if (nh > 1):
        start=ns+np+nd+nf+ng
        final=ns+np+nd+nf+ng+nh
        x[start:final] = numpy.sort(x[start:final])
    

if __name__ == "__main__":
    ns, np, nd, nf, ng, nh = get_arguments_from_cmdl(sys.argv[1:])
    nparams = ns + np + nd + nf + ng + nh

    # Load training data file
    X = numpy.loadtxt("training.dat", usecols=range(nparams))
    Y = numpy.loadtxt("training.dat", usecols=(nparams,))

    print("X shape = "+str(X.shape[0])+" x "+str(X.shape[1]))
    print("Y shape = "+str(Y.shape[0])+" x 1")
    
    # Sort parameters
    for i in range(X.shape[0]):
        sort_parameters(X[i,:], ns, np, nd, nf, ng, nh)

    # Make new training data
    data = numpy.column_stack((X,Y))

    # Check for linear dependencies
    check_for_lindep(data, ns, np, nd, nf, ng, nh)
    
    # Save new training data
    fmtstr = ""
    for v in range(nparams):
        fmtstr = fmtstr+" %10.5f"
    fmtstr = fmtstr+" %10.8f"
    numpy.savetxt("training.dat_new", data, fmt=fmtstr)
            
