import numpy as np
import subprocess
import multiprocessing as mp
import edit_plumed as ep
import matplotlib.pyplot as plt
import pandas as pd


def calc_stats(fname,op_range,boundsA,boundsB):
    
    with open(fname) as f:
        header = f.readline().split()
        assert(header[1]=="FIELDS")
        header = header[2:]
    #fh = open(fname)
    data = pd.read_csv(fname,sep='\s+',names=header,comment="#",index_col=False,skiprows=1)
    ft_idx = [i for i in range(len(data.columns)) if ("feature" in data.columns[i]) ]
    mean = data.mean()[ft_idx]
    std = data.std()[ft_idx]
    # plotd scatter matrix
    #pd.plotting.scatter_matrix(data,figsize=(12,12))
    #plt.savefig("scatter_%s.png" %  fname)
    #plt.close()


    # do histogram
    op_data = data["OP"]
    bins = np.linspace(op_range[0],op_range[1],100)
    assert(op_range[0] <= np.min(op_data))
    assert(op_range[1] >= np.max(op_data))
    hh, ee = np.histogram(op_data,bins=bins)
    hh = hh + 1.E-50
    hh = hh/np.sum(hh)
    flat = np.ones(len(hh))/len(hh)
    kld = -np.sum(flat*np.log(hh/flat))
    # calculate tranistions
    ll = len(op_data)

    fig, (ax1, ax2) = plt.subplots(2,1)
    span_time = 10 # in ps
    dt = data["time"][1]-data["time"][0]
    span = int(span_time/dt)
    rolling = np.array(op_data.rolling(span).mean())

    ax1.plot(data["time"]/1000.,op_data,lw=0.05,c='gray')
    ax1.plot(data["time"]/1000.,rolling,lw=1.25,c='k')
    
    current_state = 0
    new_state = 0
    nt = 0
    for j in range(ll):
        
        new_state = current_state
        for b in boundsA:
            if(rolling[j] >= b[0] and rolling[j] <b[1]):
                new_state = 2
        for b in boundsB:
            if(rolling[j] >= b[0] and rolling[j] <b[1]):
                new_state = 1
        if(new_state != current_state and current_state != 0):
            nt += 1
        current_state = new_state
            
            
    ax1.set_xlabel("time (ns)")
    ax1.set_ylabel("OP")

    print("%s KLD %8.4f, transitions: %d" % (fname,kld,nt))
    ax2.plot(0.5*(ee[1:]+ee[:-1]),hh,c='k')
    for b in boundsA:
        ax1.fill_between(data["time"]/1000., b[0], b[1], color='#FFB280', alpha='0.2')
        ax2.fill_betweenx([0,1.1*np.max(hh)],b[0], b[1] , color='#FFB280', alpha='0.2')
    for b in boundsB:
        ax1.fill_between(data["time"]/1000., b[0], b[1], color='#4f86f7', alpha='0.2')
        ax2.fill_betweenx([0,1.1*np.max(hh)],b[0], b[1] , color='#4f86f7', alpha='0.2')
        
    ax1.set_title("%s KLD=%5.2f NT=%d" % (fname,kld,nt) )
    
    plt.savefig("op_%s.png" % fname)
    plt.close()
    return kld,nt,mean,std



def run(idx_gen,idx_pop,weight,run_next):
    
    name = "%03d_%03d_M"% (idx_gen,idx_pop)
    pfile = ep.do_metad([op_file,features_file],name,\
                        ww=weight,mean=mean.values,std=std.values,restart=False)
        
    cmd = "%s -deffnm %s -plumed %s" % (mdrun_cmd,name,pfile)
    
    if(run_gmx):
        if(run_next<0):
            print("# RUNNING ", cmd)
            out = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
            #print("# DONE")
        else:
            oldname = "%03d_%03d_M"% (idx_gen-1,run_next)
            cmd = "cp OUTPUT_%s OUTPUT_%s " % (oldname, name)
            subprocess.check_output(cmd,shell=True)
        
    # run w static bias
    name = "%03d_%03d_S"% (idx_gen,idx_pop)
    pfile = ep.do_metad([op_file,features_file],name,\
                        ww=weight,mean=mean.values,std=std.values,restart=True)
        
    cmd = "%s -deffnm %s -plumed %s" % (mdrun_cmd,name,pfile)
    
    if(run_gmx):
        if(run_next<0):
            print("# RUNNING ", cmd)
            out = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
            #print("# DONE")
        else:
            oldname = "%03d_%03d_S"% (idx_gen-1,run_next)
            cmd = "cp OUTPUT_%s OUTPUT_%s " % (oldname, name)
            subprocess.check_output(cmd,shell=True)
                
    kld_s,nt_s,mean_s,std_s  = calc_stats("OUTPUT_%s" % name,op_range,bounds_A,bounds_B)
    return (kld_s,nt_s)



if __name__ == '__main__':
    np.random.seed(123)
    #run_gmx = False
    run_gmx = True
    generations = 100
    pop_size = 10
    elite_size = 2
    child_size = 5
    mutation_size = pop_size-2*elite_size-child_size

    # define order parameter and range

    # here, defined in plumed input file.
    # But it can be any type of classifier.
    op_file = "order_parameter.dat"
    op_range = [0,2]
    bounds_A = [[0,0.4]]
    bounds_B = [[0.9,2]]
    # define features used to construct CV
    features_file = "features.dat"

    # define other features to monitor
    # not for now

    #### OK, now we start. ###############

    # if it is not possible to bias OP run plain MD
    nsteps = 500000 
    mdrun_cmd = "gmx_mpi mdrun -s topol.tpr -nsteps %d" % (nsteps)
    all_weights = []
    pfile = ep.do_colvar([op_file,features_file],"reference")

    cmd = "%s -deffnm %s -plumed %s" % (mdrun_cmd,"reference",pfile)
    print("# RUNNING ", cmd)
    if(run_gmx):
        out = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
    print("# DONE")
    kld,nt,mean,std  = calc_stats("OUTPUT_%s" % "reference",op_range,bounds_A,bounds_B)

    nf = mean.shape[0]
    
    run_next =  [-1]*pop_size



    weights = (np.random.random((pop_size,nf))-0.5)
    pool = mp.Pool(processes=pop_size)
    fhw = open("cccc.dat","w")
    
    for idx_gen in range(generations):

        results = [pool.starmap(run,zip([idx_gen]*pop_size,range(pop_size),weights,run_next))]
        kld_list = [f[0] for f in results[0]]
        nt_list = [f[1] for f in results[0]]
        
        # create new population, elite
        idx_elite = np.argsort(kld_list)
        run_next = []
        
        new_weights = []
        print("# BEST PERFORMERS")


        for i in range(pop_size):
            if(i<elite_size):
                new_weights.append(weights[idx_elite[i]])
                run_next.append(idx_elite[i])
            print("# index %d kld %8.4f %d" % (idx_elite[i],kld_list[idx_elite[i]],nt_list[idx_elite[i]]))
            ss1 = " %d kld %8.4f %d " % (idx_elite[i],kld_list[idx_elite[i]],nt_list[idx_elite[i]])
            ss1 += " ".join(["%8.4e" % r for r in weights[idx_elite[i]] ])
            fhw.write(ss1+"\n")
            fhw.flush()
        # create new population, crossover
        for i in range(child_size):

            # pick two random from best half
            aa = idx_elite[:int(0.5*len(kld_list))]
            np.random.shuffle(aa)
            parent_1 = int(aa[0])
            parent_2 = int(aa[1])

            # select random mutation
            mut_point = np.random.randint(nf)
            child = list(weights[parent_1][:mut_point]) + list(weights[parent_2][mut_point:])
            new_weights.append(child)
            run_next.append(-1)

        # create new_population, single point mutation
        for i in range(elite_size):
            # pick elite 
            idx_parent = idx_elite[np.random.randint(elite_size)]
            # pick random position
            idx_mut = np.random.randint(nf)
            
            # pick random weight
            new_weight = weights[idx_parent]
            new_weight[idx_mut] =  (np.random.random()-0.5)

            new_weights.append(new_weight)
            run_next.append(-1)
        
        # create new_population, random mutation
        for i in range(mutation_size):
            # pick random parent 
            #idx_parent = np.random.randint(pop_size)
            # pick random position
            #idx_mut = np.random.randint(nf)
            # pick random weight
            new_weight =  (np.random.random(nf)-0.5)
            #weights[idx_parent]
            #new_weight[idx_mut] = 

            new_weights.append(new_weight)
            run_next.append(-1)
            

        weights = np.copy(new_weights)
        #            print()
        
        all_weights.append(weights)
        
        #print(" NEXT iteration")
        #for k in range(len(weights)):
        #    print(k,weights[k])
            
