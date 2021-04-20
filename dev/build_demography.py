import msprime
import ecophylo
import numpy as np

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
        [0.6,0]]]
migr = [1, 0.5, 0.8] 
changemigr = [0, 100, 200]

rates = [[0,-0.002], [0,0]]
changerates = [[0,100], [600, 80]]

#def build_demography(samples, comsizes, changetimes, rates = None, changerates = None, migr = None, changemigr = None):
if True :
    """
    """
    
    # IDIOT PROOF HERE
    # check that there are as many changetimes[pop] as comsize[pop]
     

    npop = len(samples)
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
        if len(changerates[pop]) > 1:
            for i in range(len(changerates[pop][1:])):
                demography.add_population_parameters_change(time = changerates[pop][i+1] , growth_rate=rates[pop][i+1], population= pop_ids[pop])

    # initialize migration matrix    
    if migr is not None:
        m = np.array(migr)
        dim = m.shape
        if np.sum(m) == 0 :
            print("ERROR: migration matrices cannot all be empty")

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
                            if migr[t+1][row][col] == 0 :
                                continue
                            demography.add_migration_rate_change(time = changemigr[t], rate = migr[t+1][row][col], source=pop_ids[row], dest=pop_ids[col])



demography.sort_events()
print(demography.debug())
#ts = msprime.sim_ancestry(demography = demography, samples = sampleset)
#tree = ts.first()
#print(tree.draw(format = "unicode"))
#range(len(changetimes[pop][1:]))