# -*- coding: utf-8 -*-

# def foo (x= None, **kwargs):
#     for count, thing in kwargs.items():
#         # ["mot_cle"]
#         print('{0}. {1}'.format(count, thing))
#     return x

# def foo (x= None, *args):
#     for count, thing in enumerate(args):
#         # ["mot_cle"]
#         print('{0}. {1}'.format(count, thing))
#     return x


# prior_locate = [["com_size", 0, 1, [50, 150, "unif"]],
#                 ["com_size", 1, 1, [300, 500, "unif"]]]    

import pandas as pd
import numpy as np
from ecophylo import check_params
from ecophylo import simulate
from ecophylo import getAbund
from ecophylo import getDeme

# t = simulate(samples = [2, 2], com_size = [[20, 50, 100], [20, 50]], mu = 1, 
# init_rates = None, 
# changetime = [[0, 20, 30], [0, 20]] , mrca = None, 
# migr = None, migr_time = None, vic_events = None,
# verbose = True, seed = 42)
# TODO : expect error here !!

def dosimuls(nsim, samples, com_size, mu, init_rates = None, changetime = None,
             mrca = None, migr = 1, migr_time = None,
             vic_events = None, 
             verbose = False,
             output = ['Params'], # Params, Sumstat, Tree
             file_name = None, seed = None):

    # Idiotproof dosimul parameters
    # nsim
    if not isinstance(nsim, (float,int)) or nsim < 1 :
        raise ValueError("nsim must be a int value above 0")
    if isinstance(nsim, float) : nsim = int(nsim)
    # output
    if not isinstance(output, list) or any([not isinstance(x, str) for x in output]) :
        raise ValueError("output must be a list of string")
    if len(output) > 3 or any([ x not in ["Params", "Sumstat", "Trees"] for x in output]):
        raise ValueError("output must contain only 'Params', 'Sumstat' or/and 'Trees'")
    # file_name
    if file_name is not None and not isinstance(file_name, str) :
        raise ValueError("file_name must be string value.")
            

    # register None values before they are initiated
    Nonedef = [init_rates, changetime, mrca, migr_time, vic_events]
    for i in range(len(Nonedef)):
        if Nonedef[i] is not None : Nonedef[i] = True
        else : Nonedef[i] = False
    
    # IDIOTPROOF and Register locations
    samples, com_size, mu, init_rates, changetime, mrca, migr, migr_time, \
                    vic_events, verbose, seed, prior_locate = check_params(
        samples = samples, com_size = com_size, mu = mu, init_rates = init_rates, 
        changetime = changetime, mrca = mrca, migr = migr, 
        migr_time = migr_time, vic_events = vic_events, 
        verbose = verbose, seed = seed, prior_locate = 'naive'
    )
    npop = len(samples)
    priors = [prior_locate[p][0] for p in range(len(prior_locate))]
    print("prior_locate :",prior_locate) # TODO : remove

    ##########################################################################
    ####                    SETUP DATAFRAME                               ####
    ##########################################################################
    # Creating parameter dataframe
    # if prior locate is empty deal with this later
    comments = ("Simulating eco-evolutionary dynamics over the following" +
               " parameters ranges:")
    params = pd.DataFrame()

    prior_names = []
    for prior in prior_locate :
        if prior[0] == "samples": prior_names.append(f'samples_pop{prior[1]}')
        if prior[0] == "mu": prior_names.append('mu')
        if prior[0] == "com_size":
            prior_names.append(f'com_size_pop{prior[1]}_t{prior[2]}')
        if prior[0] == "changetime":
            prior_names.append(f'time_pop{prior[1]}_t{prior[2]}')
        if prior[0] == "init_rates":
            prior_names.append(f'rate_pop{prior[1]}_t{prior[2]}')
        if prior[0] == "migr": prior_names.append(f'migr_t{migr_time[prior[1]]}')
    ## samples
    col_samples = ['samples_pop{}'.format(i) for i in range(len(samples))]
    for j in range(len(samples)):
        if col_samples[j] in prior_names: params[col_samples[j]] =  [None]*nsim
        else: params[col_samples[j]] = [samples[j]]*nsim
    ## mu
    if "mu" in priors: params["mu"] = [None]*nsim 
    else: params["mu"] = [mu]*nsim
    ## com_size
    for j in range(len(com_size)):
        for jj in range(len(com_size[j])):
            tmp_name = f'com_size_pop{j}_t{jj}'
            if tmp_name in prior_names: params[tmp_name] = [None] * nsim
            else : params[tmp_name] = [com_size[j][jj]] * nsim
    ## changetime
    if Nonedef[1]:
        for i in range(npop):
            for ii in range(len(changetime[i])):
                if ii != 0:
                    tmp_name = f'time_pop{i}_t{ii}'
                    if tmp_name in prior_names: params[tmp_name] = [None] * nsim
                    else : params[tmp_name] = [changetime[i][ii]] * nsim
    ## init_rates
    if Nonedef[0]:
        for i in range(npop):
            for ii in range(len(init_rates[i])):
                tmp_name = f'rates_pop{i}_t{ii}'
                if tmp_name in prior_names: params[tmp_name] = [None] * nsim
                else : params[tmp_name] = [init_rates[i][ii]] * nsim
    ## migr and migr_time
    if migr is not None and isinstance(migr[0], (int, float)):
        for i in range(len(migr)):
            tmp_name = f'migr_t{migr_time[i]}'
            if tmp_name in prior_names: params[tmp_name] =  [None]*nsim
            else: params[tmp_name] = [migr[i]]*nsim
        # # This should work for matrices but we don't want to fuck our brains
        # col_migr = []
        # migr_val = []
        # for mt in range(len(migr)):
        #     if isinstance(migr[mt], (int, float)):
        #         migr[mt] = np.ones((npop,npop))*migr[mt]
        #         np.fill_diagonal(migr[mt], 0)
        #     x = np.array(migr[mt])
        #     XX,YY = np.meshgrid(np.arange(x.shape[1]),np.arange(x.shape[0]))
        #     table = np.vstack((x.ravel(),XX.ravel(),YY.ravel())).T
        #     col_migr.extend([f'migr_pop{int(t[1])}_pop{int(t[2])}_t{mt}' for t in table if np.unique(t[1:]).size != 1])
        #     migr_val.extend([t[0] for t in table if np.unique(t[1:]).size != 1])
        # for j in range(len(col_migr)):
        #     if 'migr' in priors:
        #         params[col_migr[j]] =  [None]*nsim
        #     else: 
        #         params[col_migr[j]] = [migr_val[j]]*nsim

    ## vic_events
    # if Nonedef[4]:
    #     print("you are fucked") # TODO : maybe remove this
    # print(params)

    result = list()
    abund = list()
    diver = list()
    trees = ""
    print("\n****SIMULATING***************\n\n")
    failed = 0
    i = 0
    while i < nsim and failed < 100:
        # SIMULATE
        phylo = simulate(
            samples = samples, com_size = com_size, mu = mu, 
            init_rates = init_rates, changetime = changetime, 
            mrca = mrca, migr = migr, migr_time = migr_time, 
            vic_events = vic_events, verbose = False, seed = seed, force = True
        )
        #################################################################
        ####                    CHECK TREE                           ####
        #################################################################
        # check tree
        sp  = phylo.get_leaf_names()
        if len(sp) > 2 and len(set(sp)) > 1 :
            for ii in range(len(prior_locate)):
                tmp_p = prior_locate[ii]
                if tmp_p[0] == "samples":
                    print('not done') # TODO : samples 
                if tmp_p[0] == "mu":
                    params.loc[i,'mu'] = mu
                if tmp_p[0] == "com_size":
                    params.loc[i,(f'com_size_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        com_size[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == "changetime":
                    params.loc[i,(f'time_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        changetime[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == "init_rates":
                    params.loc[i,(f'rates_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        init_rates[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == 'migr':
                    params.loc[i,(f'migr_t{migr_time[tmp_p[1]]}')] = \
                        migr[tmp_p[1]]
            # Save sumstat
            abund.append(getAbund(phylo))
            diver.append(getDeme(phylo, div = True))
            # Save tree
            trees += phylo.write() + "\n"
            failed = 0
            i+=1
        else :
            failed += 1      
        #################################################################
        ####                    RESAMPLING                           ####
        #################################################################
        samples, com_size, mu, init_rates, changetime, mrca, migr, migr_time, \
                    vic_events, verbose, seed, prior_locate = check_params(
            samples = samples, com_size = com_size, mu = mu, init_rates = init_rates, 
            changetime = changetime, mrca = mrca, migr = migr, 
            migr_time = migr_time, vic_events = vic_events, 
            verbose = verbose, seed = seed, prior_locate = prior_locate
        )

    if failed >= 100 : raise ValueError("Too many simulations have failed")

    if verbose:
        print(comments)

    # abund
    abund_pd = np.zeros((nsim, max(len(x) for x in abund)))
    for k, j in enumerate(abund):
        abund_pd[k][0:len(j)] = j
    abund_pd = pd.DataFrame(abund_pd)
    # diver
    diver_pd = pd.DataFrame(diver)
    print(['alpha'+ str(x) for x in range(diver_pd.shape[1] )])
    diver_pd.columns = ['alpha'+ str(x) for x in range(diver_pd.shape[1] )]
    print(diver_pd)
    sumstat = pd.concat([abund_pd.reset_index(drop=True), diver_pd], axis=1)

    print(sumstat)

    if file_name is not None :
        saved = ""
        if "Params" in output:
            saved += params.to_string()
        if "Sumstat" in output :
            saved += "\n###\n" + sumstat.to_string()
        if "Trees" in output : 
            saved += "\n###\n" + trees

        print('\nSimulation saved in :', file_name, '\n')
        f = open(file_name, "w")
        f.write(saved)
        f.close()

    result = [params, sumstat, trees]
    return result

print("\n\n* test 1\n")
print(dosimuls(nsim = 5, samples = [10],
    com_size = [[500]],
    mu = 0.1, # will fail if 0.0001
    seed = 42))
# print("\n\n* test mu\n")
# print(dosimuls(nsim = 5, samples = [10],
#     com_size = [[500]],
#     mu = [0.001, 0.005, 'uniform']))
print("\n\n* test 2\n")
print(dosimuls(nsim = 5, samples = [10],
    com_size = [[[1500, 5000, "uniform"], [1500, 5000, "uniform"]]],
    mu = 0.001, changetime = [0, 500], seed = None,  file_name = 'tmp.txt'))
print("\n\n* test 3\n")
print(dosimuls(nsim = 5, samples = [10, 9],
    com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
    mu = 0.001, migr = 0.5,
    changetime = [[0, 200], [0, 100, 500]],
    seed = 42))
print("\n\n* test migr\n")
print(dosimuls(nsim = 5, samples = [10, 9],
    com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
    mu = 0.001, migr = [[0.2,0.5,"uniform"],[0.2,0.5,"uniform"]], 
    file_name = 'tmpmigr.txt', migr_time = [0, 200],
    output = ["Params", "Sumstat", "Trees"],
    changetime = [[0, 200], [0, 100, 500]]))

# print("\n\n* test migr\n")
# print(dosimuls(nsim = 5, samples = [10, 9],
#     com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = [[[0, 0.5], [0.2, 0]]],
#     changetime = [[0, 200], [0, 100, 500]], file_name = 'tmp.txt',
#     seed = 42))


# print(dosimuls(nsim = 5, samples = [10, 10],
#     com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = 0.5,
#     changetime = [[0, 200], [0, 100, 500]]))

# print(check_params(
#     samples = [10, 10],
#     com_size = [[500, 1000], [200, [150, 500, "uniform"], 600]],
#     mu = 0.001,
#     changetime = [[0, 200], [0, 100, 500]],
#     seed = 42,
#     prior_locate = "naive"))

# print(check_params(
#     samples = [10, 10],
#     com_size = [[500, 1000], [200, 373, 600]],
#     mu = 0.001,
#     changetime = [[0, 200], [0, 100, 500]],
#     seed = 42,
#     prior_locate = [['com_size', 1, 1, [150, 500, 'uniform']]]))

# print(check_params(
#     samples = [10, 10],
#     com_size = [[500, 1000], [200, [150, 500, "uniform"], 600]],
#     mu = 0.001,
#     changetime = [[0, 200], [0, 100, 500]],
#     seed = 42,
#     prior_locate = None))
