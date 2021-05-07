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
from ecophylo import check_params
from ecophylo import simulate
from ecophylo import getAbund
from ecophylo import getDeme

def dosimuls(nsim, samples, com_size, mu, init_rates = None, changetime = None,
             mrca = None, migr = 1, migr_time = None,
             vic_events = None, 
             verbose = False, sumstat = None, result = ['Params'], 
             file_name = None, seed = None):

    # Idiotproof dosimul parameters
    # sumstat
    # result
    # file_name
    if file_name is not None and not isinstance(file_name, str):
        raise ValueError("file_name must be string value.")


    # register None value before they are initiated
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
    prior_names = [prior_locate[p][0] for p in range(len(prior_locate))]
    print("prior_locate :",prior_locate) # TODO : remove

    ##########################################################################
    ####                    SET DATAFRAME                                 ####
    ##########################################################################
    # Creating parameter dataframe
    # if prior locate is empty deal with this later
    comments = ("Simulating eco-evolutionary dynamics over the following" +
               " parameters ranges:")
    df = pd.DataFrame()
    params = {"samples" : samples, "com_size" : com_size, "mu" : mu,
              "init_rates": init_rates, "changetime" : changetime,
              "migr" : migr, "migr_time": migr_time, "vic_events": vic_events}

    params = {k: v for k, v in params.items() if v is not None}
    
    # cst_params = {k: v for k, v in params.items() if k not in prior_names} # useless
    prior_names = [f'{prior_names[p]}_pop{prior_locate[p][1]}_t{prior_locate[p][2]}' for p in range(len(prior_locate))]

    ## samples
    col_samples = ['samples_pop{}'.format(i) for i in range(len(samples))]
    for j in range(len(samples)):
        if col_samples[j] in prior_names: df[col_samples[j]] =  [None]*nsim
        else: df[col_samples[j]] = [samples[j]]*nsim

    ## com_size
    col_com_size = [f'com_size_pop{index1}_t{index2}' for index1,value1 in enumerate(com_size) for index2,value2 in enumerate(value1)]
    lst_col_com_size = [col_com_size[i:i + len(com_size)] for i in range(0, len(col_com_size), len(com_size))]
    # print(lst_col_com_size)
    for j in range(len(lst_col_com_size)):
        for jj in range(len(lst_col_com_size[j])):
            if lst_col_com_size[j][jj] in prior_names:
                df[lst_col_com_size[j][jj]] =  [None]*nsim
            else:
                df[lst_col_com_size[j][jj]] =  [com_size[j][jj]]*nsim

    ## mu
    df["mu"] = [mu]*nsim

    ## init_rates
    if Nonedef[0]:
        col_init_rates = [f'rate_pop{index1}_t{index2}' for index1,value1 in enumerate(com_size) for index2,value2 in enumerate(value1)]
        for col in col_init_rates:
            df[col] = ""
           
    ## changetime
    if Nonedef[1]:
        col_times = [f'time_pop{index1}_t{index2}' for index1,value1 in enumerate(com_size) for index2,value2 in enumerate(value1)]
        for col in col_times:
            df[col] = ""

    # ## migr
    # if Nonedef[]:
    #     print("find a way to do this") # TODO : maybe remove this

    # ## migr_time
    # if all(Nonedef[]):
    #     col_migrtime = ['migr_t{}'.format(i) for i in range(len(migr))]
    #     for col in col_migrtime:
    #         df[col] = ""
           
    ## vic_events
    if Nonedef[4]:
        print("you are fucked") # TODO : maybe remove this
  
    print(df)
    ##########################################################################
    ####                    SIMULATING                                    ####
    ##########################################################################

    print("\n SIMULATING \n\n============\n")
    failed = 1
    i = 0
    while i < nsim and failed <= 100:

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
        if True : # TODO : if tree ok
            # add parameters
            for ii in range(len(prior_locate)):
                tmp_p = prior_locate[ii]
                if tmp_p[0] == "samples":
                    print('not done')
                if tmp_p[0] == "com_size":
                    df.loc[i,(f'com_size_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        com_size[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == "mu":
                    print('not done')
                if tmp_p[0] == "init_rates":
                    print('not done')
                if tmp_p[0] == "changetime":
                    print('not done')
                if tmp_p[0] == "migr":
                    print('not done')

            # Sumstat

            failed = 1
            i+=1
        else :
            failed += 1
            # maybe resample value ?
      
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

    if failed > 100 :
        raise ValueError("Too many simulations have failed")

    if file_name is not None :
        saved = df.to_string()

        print('need to save the file')
        f = open(file_name, "w")
        f.write(saved)
        f.close()

    return df

print("\n\n* test 1\n")
print(dosimuls(nsim = 5, samples = [10],
    com_size = [[500]],
    mu = 0.001, 
    seed = 42))
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


