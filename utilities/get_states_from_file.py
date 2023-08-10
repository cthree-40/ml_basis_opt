import sys, getopt

def get_arguments_from_cmdl(argv):
    nstates = 0
    opts, args = getopt.getopt(argv,"n:m:")
    for opt, arg in opts:
        if (opt == "-n"):
            nstates = arg
        elif (opt == "-m"):
            molec = arg
        else:
            print ("Unknown option")
            sys.exit()
    return (int(nstates), molec)

def create_states_file(molec, nstates):
    qfile = open(molec+".output", "r")
    qf_lines = (qfile.read().splitlines())
    qfile.close()
    nlines = len(qf_lines)
    estart_line = 0
    for i in range(nlines):
        if "==== Final CI Energy ====" in str(qf_lines[i]):
            estart_line = i
            break
    estart_line = estart_line + 1
    energy = []
    for i in range(nstates):
        e_line = qf_lines[estart_line + i]
        e_line = e_line.replace(" CI Energy (au) Root #   "+str(i),"").strip()
        energy.append(e_line)

    efile = open(molec+"_states.data", "w")
    for i in range(nstates):
        efile.write("%s\n" % energy[i])
    efile.close()

if __name__ == "__main__":
    nstates, molec = get_arguments_from_cmdl(sys.argv[1:])
    create_states_file(molec, nstates)
