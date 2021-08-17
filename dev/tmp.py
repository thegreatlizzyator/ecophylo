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

import ecophylo
# mu = 0.05
# samples = [3,3,3]
# deme_sizes = [[200], [200, 300, 400], [200]]
# change_time = [[0], [0, 1000, 2000], [0]]
# splits = [[[0,1], 1, 1],
#           [[1,2], 1, 2]]

# tree = ecophylo.simulate(samples = samples,
#                          deme_sizes= deme_sizes, 
#                          mu = mu, 
#                          changetimes= changetimes,
#                          splits = splits, seed= 42)

# import pandas as pd
# import numpy as np
# from ecophylo import check_params
# from ecophylo import simulate
# from ecophylo import getAbund
# from ecophylo import getDeme
# from ecophylo import dosimuls

# t = simulate(samples = [5, 5], deme_sizes = [[2e3, 4e3], [2e3]], mu = 1,
#              gr_rates = None, changetimes = [[0, 300], [0]] , mrca = None, 
#              migr = None, migr_times = None, splits = [[[0, 1], 0, 1]],
#              verbose = False, seed = 42)
# print(t)

# lab = {u: 'sp'+str(u) for u in range(len(t.get_leaf_names())) }

# t = simulate(samples = [2, 2], com_size = [[20, 50, 100], [20, 50]], mu = 1, 
# init_rates = None, 
# changetime = [[0, 20, 30], [0, 20]] , mrca = None, 
# migr = None, migr_time = None, vic_events = None,
# verbose = True, seed = 42)
# TODO : expect error here !!


# print("\n\n* test 1\n")
# print(dosimuls(nsim = 5, samples = [10],
#     deme_sizes = [[500]],
#     mu = 0.1, # will fail if 0.0001
#     seed = 42))
# print("\n\n* test mu\n")
# print(dosimuls(nsim = 5, samples = [10],
# #     com_size = [[500]],
# #     mu = [0.001, 0.005, 'uniform']))
# print("\n\n* test 2\n")
# print(dosimuls(nsim = 5, samples = [10],
#     deme_sizes = [[[1500, 5000, "uniform"], [1500, 5000, "uniform"]]],
#     mu = 0.001, changetimes = [0, 500], seed = 42))
# print("\n\n* test 3\n")
# print(dosimuls(nsim = 5, samples = [10, 9],
#     deme_sizes = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = 0.5,
#     changetimes = [[0, 200], [0, 100, 500]],
#     seed = 42))
# print("\n\n* test migr\n")
# print(dosimuls(nsim = 5, samples = [10, 9],
#     deme_sizes = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = [[0.2,0.5,"uniform"],[0.2,0.5,"uniform"]], 
#     migr_times = [0, 200],
#     output = ["Params", "Sumstat", "Trees"],
#     changetimes = [[0, 200], [0, 100, 500]]))
# print("\n\n* test migr\n")
# print(dosimuls(nsim = 5, samples = [10, 9],
#     deme_sizes = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = [[[0, 0.5], [0.2, 0]]],
#     changetimes = [[0, 200], [0, 100, 500]],
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


t = ecophylo.dosimuls(nsim = 10, samples = 3063, deme_sizes = [[6e5, 1e11]], 
mu = [1e-7, 1e-9, "uniform"], 
gr_rates = None, changetimes = [[0, 8e4]] , mrca = None, migr = 1, 
migr_times = None, file_name = "tmp.lol", verbose = True)

print(t[0].loc[0:1,"mu"])
print(t[0].loc[0,"mu"])

# import pandas as pd
# d = {'col1': [1, 2], 'col2': [6.000822964516544e-08, 4.034576430589255e-08]}
# df = pd.DataFrame(d)
# print(df)
# print(df.to_string())