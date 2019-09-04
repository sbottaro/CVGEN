

def calc_obs(fname):

    stri = ""
    labels = []
    fh = open(fname)
    
    for line in fh:
        if(len(line.split())==0): continue
        # add line
        stri += line
        
        # get labels to print
        ss = line.split(":")
        if(len(ss) == 2 and "feature" in ss[0]):
            labels.append(ss[0])
        if("LABEL=OP" in line):
            labels.append("OP")
            
    stri += "\n"

    if(len(labels)==0):
        print("# ERROR: nothing to print in %s " % fname)
        exit(1)

            
    return stri,labels



def do_colvar(files,suffix,stride=100):

    stri = ""
    labs = ""
    for f in files:
        ss, ll = calc_obs(f) 
        stri += ss
        labs = "%s,%s" % (labs,",".join(ll))

    stri += "PRINT FILE=OUTPUT_%s STRIDE=%d ARG=%s" % (suffix,stride,labs)
    fname = "plumed_%s.dat" % suffix
    fhw = open(fname,"w")
    fhw.write(stri)
    fhw.close()
    return fname


def create_cv(fname,ww,mean,sigma):
    
    for line in fh:
        if(len(line.split())==0): continue
        # add line
        stri += line
        
        # get labels to print
        ss = line.split(":")
        if(len(ss) == 2 and "feature" in ss[0]):
            labels.append(ss[0])
    assert(len(ww)==len(sigma))
    assert(len(ww)==len(mean))
    
    stri = "COMBINE ...\n LABEL=CV\n PERIODIC=NO\n"
    stri += "ARG=%s \n" % (",".join(labels))
    #stri += "COEFFICIENTS=%s \n" % (",".join(["%8.4e" % ww[j]/sigma[j] for j  in range(len(ww))]))
    #stri += "PARAMETERS=%s \n" % (",".join(["%8.4e" % mean[j] for j  in range(len(ww))]))
    stri += "COEFFICIENTS=%s \n" % (",".join(["%8.4e" % ww[j] for j  in range(len(ww))]))
    #stri += "PARAMETERS=%s \n" % (",".join(["%8.4e" % mean[j] for j  in range(len(ww))]))

    stri += "... COMBINE"


def do_metad(files,suffix,ww,mean,std,restart,stride=100):

    stri = ""
    if(restart):
        stri += "RESTART\n"

    labels = []
    labs = ""
    for f in files:
        ss, ll = calc_obs(f) 
        stri += ss
        labs = "%s,%s" % (labs,",".join(ll))
        for el in ll:
            if("feature" in el): labels.append(el)


    cc = ww/std

    stri += "COMBINE ...\n LABEL=CV\n PERIODIC=NO\n"
    stri += "ARG=%s \n" % (",".join(labels))
    stri += "COEFFICIENTS=%s \n" % (",".join(["%8.4e" % cc[j] for j  in range(len(ww))]))
    stri += "PARAMETERS=%s \n" % (",".join(["%8.4e" % mean[j] for j  in range(len(ww))]))
    stri += "... COMBINE\n"
    
    labs+= ",CV"
    
    stri += "PRINT FILE=OUTPUT_%s STRIDE=%d ARG=%s RESTART=NO\n" % (suffix,stride,labs[1:])
    if(restart):
        stri += "METAD ARG=CV PACE=500000000 HEIGHT=0.0 SIGMA=0.1 FILE=bias_%s\n" % (suffix.replace("S","M"))
    else:
        stri += "METAD ARG=CV PACE=100 HEIGHT=0.5 SIGMA=0.1 FILE=bias_%s\n" % (suffix)

    fname = "plumed_%s.dat" % suffix
    fhw = open(fname,"w")
    fhw.write(stri)
    fhw.close()
    return fname

