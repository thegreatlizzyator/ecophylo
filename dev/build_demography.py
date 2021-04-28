import msprime
import ecophylo
import numpy as np
import sys 

samples= [5, 5]
comsizes = [[10,20],[40,50,60]]
changetimes = [[0,100], [0,300,400]]

migr = [[[0,0.1],
        [0.2,0]]
        ,
        [[0,0.3],
        [0.4,0]]
        ,
        [[0,0.5],
        [0.6,0]]
        ,
        [[0,0.7],
        [0.6,0]]]

migr = [1, 0.5, 0.8, 0.7]

changemigr = [0, 100, 200, 300]

rates = [[0,-0.002], [0,0]]

#def build_demography(samples, comsizes, changetimes, rates = None, changerates = None, migr = None, changemigr = None):
if True :
    """
    """
    npop = len(samples)

    #sample idiot proof
    isint_samp = [isinstance(s, int) for s in samples]
    ispos_samp = [s>0 for s in samples]
    if not all(isint_samp):
        sys.exit("sample sizes should all be ints")
    if not all(ispos_samp):
        sys.exit("sample sizes should all be positive")
    
    # past sizes idiot proof
    if changetimes is not None:
        if len(comsizes) != npop:
            sys.exit("there should be as many elements in comsizes as there are demes")
        for p in range(npop):
            if len(comsizes[p]) != len(changetimes[p]):
                sys.exit("there should be as many past sizes as there are epochs in changetime")
        for sizes in comsizes:
                isint_sizes = [isinstance(s, int) for s in sizes]
                ispos_sizes = [s>0 for s in sizes]
                if not all(isint_sizes) :
                    sys.exit("all past sizes should be ints")
                if not all(ispos_sizes):
                    sys.exit("all past sizes should be strictly positive")

    # growth rates idiot proof
    if rates is not None: 
        if len(rates) != npop:
            sys.exit("there should be as many elements in rates as there are demes")
        for p in range(npop):
            if len(rates[p]) != len(changetimes[p]):
                sys.exit("there should be as many past growth rates as there are epochs in changetime")
        for gr in rates:
                isint_rates = [isinstance(r, (int, float)) for r in gr]
                if not all(isint_sizes) :
                    sys.exit("all past growth rates should be ints or floats")
        

    # migration idiot proof
    if migr is not None:
        if npop == 1 :
            sys.warn("no migration matrix is needed for a single deme")
            migr = None
        m = np.array(migr)
        dim = m.shape

        if np.sum(m) == 0 :
            sys.exit("migration matrices cannot all be empty")

        if dim[0] != len(changemigr):
            sys.exit("there should be as many symmetric migration rates as there are times in changemigr")
        if len(dim) > 1:
            if dim[1] != dim[2] or dim[1] != npop or dim[2] != npop:
                sys.exit("custom migration matrices should be of size ndeme x ndeme")
            for mat in migr:
                isint_migr = [isinstance(i, (float,int)) for i in mat]
                ispos_migr = [i>=0 or i<= npop for i in mat]
                if not all(isint_migr) :
                    sys.exit("found custom migration matrix that is not made of ints or floats")
                if not all(ispos_migr):
                    sys.exit("found custom migration matrix with negative migration rates or greater than the number of demes")
        else:
            isint_migr = [isinstance(i, (float,int)) for i in migr]
            ispos_migr = [i>=0 or i<= npop for i in migr]
            if not all(isint_migr):
                sys.error("migration rates should be either ints or floats")
            if not all(ispos_migr): 
                sys.error("migration rate should be positive (or zero) and not exceed the number of demes")
    

    pop_ids = ["pop_" + str(p) for p in range(npop)]
    demography = msprime.Demography()

    # build SampleSet
    sampleset = [msprime.SampleSet(samp, pop, ploidy=1) for samp,pop in zip(samples, pop_ids)]
    
    # initialize population configurations
    for pop in range(npop):
        demography.add_population(initial_size= comsizes[pop][0], growth_rate=rates[pop][0])
        
        # if population sizes have fluctuated in the past:
        if len(changetimes[pop]) > 1:
            for i in range(len(changetimes[pop][1:])):
                demography.add_population_parameters_change(time = changetimes[pop][i+1] , initial_size=comsizes[pop][i+1], population= pop_ids[pop])
        
        # if population growth rates have fluctuated in the past:
        if len(changetimes[pop]) > 1:
            for i in range(len(changetimes[pop][1:])):
                demography.add_population_parameters_change(time = changetimes[pop][i+1] , growth_rate=rates[pop][i+1], population= pop_ids[pop])

    # initialize migration matrix    
    if migr is not None:
        if len(dim) == 1 :
            # symmetric migration matrix
            demography.set_symmetric_migration_rate(populations = range(npop), rate = migr[0])

            # if symmatric migration rate has changed in the past:
            if len(migr) > 1:
                for m in range(len(migr[1:])):
                    demography.add_migration_rate_change(time = changemigr[m+1], rate = migr[m+1])

        if len(dim) > 2 :
            # custom migration matrix
            for row in range(npop):
                for col in range(npop):
                    if migr[0][row][col] == 0 :
                        continue
                    demography.set_migration_rate(source = pop_ids[row], dest = pop_ids[col], rate = migr[0][row][col])

            # if migration matrix has changed in the past
            if len(migr)>1:
                for t in range(len(changemigr)):
                    for row in range(npop):
                        for col in range(npop):
                            if migr[t][row][col] == 0 :
                                continue
                            demography.add_migration_rate_change(time = changemigr[t], rate = migr[t][row][col], source=pop_ids[row], dest=pop_ids[col])



demography.sort_events()
#print(demography.debug())
#ts = msprime.sim_ancestry(demography = demography, samples = sampleset)
#tree = ts.first()
#print(tree.draw(format = "unicode"))
#range(len(changetimes[pop][1:]))