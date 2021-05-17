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

t = simulate(samples = [5, 5], com_size = [[2e3, 4e3], [2e3]], mu = 1,
             init_rates = None, changetime = [[0, 300], [0]] , mrca = None, 
             migr = None, migr_time = None, vic_events = [[[0, 1], 0, 1]],
             verbose = False, seed = 42)
print(t)

lab = {u: 'sp'+str(u) for u in range(len(t.get_leaf_names())) }

# t = simulate(samples = [2, 2], com_size = [[20, 50, 100], [20, 50]], mu = 1, 
# init_rates = None, 
# changetime = [[0, 20, 30], [0, 20]] , mrca = None, 
# migr = None, migr_time = None, vic_events = None,
# verbose = True, seed = 42)
# TODO : expect error here !!

# print("\n\n* test 1\n")
# print(dosimuls(nsim = 5, samples = [10],
#     com_size = [[500]],
#     mu = 0.1, # will fail if 0.0001
#     seed = 42))
# # print("\n\n* test mu\n")
# # print(dosimuls(nsim = 5, samples = [10],
# #     com_size = [[500]],
# #     mu = [0.001, 0.005, 'uniform']))
# print("\n\n* test 2\n")
# print(dosimuls(nsim = 5, samples = [10],
#     com_size = [[[1500, 5000, "uniform"], [1500, 5000, "uniform"]]],
#     mu = 0.001, changetime = [0, 500], seed = None,  file_name = 'tmp.txt'))
# print("\n\n* test 3\n")
# print(dosimuls(nsim = 5, samples = [10, 9],
#     com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = 0.5,
#     changetime = [[0, 200], [0, 100, 500]],
#     seed = 42))
# print("\n\n* test migr\n")
# print(dosimuls(nsim = 5, samples = [10, 9],
#     com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]],
#     mu = 0.001, migr = [[0.2,0.5,"uniform"],[0.2,0.5,"uniform"]], 
#     file_name = 'tmpmigr.txt', migr_time = [0, 200],
#     output = ["Params", "Sumstat", "Trees"],
#     changetime = [[0, 200], [0, 100, 500]]))

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
