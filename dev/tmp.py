# -*- coding: utf-8 -*-

# import numpy as np
# migr = [[[0,0.1],
#         [0.2,0]]
#         ,
#         [[0,0.3],
#         [0.4,0]]
#         ,
#         [[0,0.3],
#         [0.4,0]]
#         ,
#         [[0,0.3],
#         [0.4,0]]
#         ,
#         [[0,0.5],
#         [0.6,0]]]



# -*- coding: utf-8 -*-
# """
# Created on Wed Apr 21 14:36:39 2021

# @author: barthele
# """
# import msprime
# import collections

# samples = [3,3,3]
# com_size = [[200], [200, 300, 400], [200]]
# changetime = [[0], [0, 1000, 2000], [0]]
# vic_events = [[1000, [0,1], 1],
#               [2000, [1,2], 1]]

# catch = [[1000, [0,1], 1],
#          [2000, [1,2], 2]] #TODO: try and catch this

# def flatten(x):
#     if isinstance(x, collections.Iterable):
#         return [a for i in x for a in flatten(i)]
#     else:
#         return [x]

#  #def vic(samples, com_size, changetime, vic_events):
# if True:
#     npop = len(samples)
#     demography_test = msprime.Demography()
#     #pop_ids = ["pop_" + str(p) for p in range(npop)]
    
#     # initialize population configurations
#     for pop in range(npop):
#         demography_test.add_population(initial_size= com_size[pop][0])

#     if vic_events is not None:
#         # idiot proof 
#         if any([len(v)!=3 for v in vic_events]) :
#             raise ValueError('all elements in vic_events should be  lists of lenght 3 ')
        
#         if any([len(v[1])!=2 for v in vic_events]):
#             raise ValueError('second element of vic_events should be a list of 2 deme ids')
        
#         if not all([isinstance(v, (int)) for v in flatten(vic_events)]):
#             raise ValueError("all elements of vic_events should be ints")
        
#         if any([t<0 for t in [v[0] for v in vic_events]]):
#             raise ValueError("all times in _vic_events should be strictly positive")
        
#         if any([test not in flatten(changetime) for test in [v[0] for v in vic_events]]):
#             raise ValueError("split times in vic_events should also appear in changetime")

#         if set([p for p in range(len(samples))]) != set(sum([flatten(v[1:]) for v in vic_events],[])):
#             raise ValueError("Split events do not match provided deme information")
        
#         vic_dates = [v[0] for v in vic_events]
#         if vic_dates != sorted(vic_dates):
#             raise ValueError("Split dates should be provided in chronological order")
        
#         if any([v[2] not in v[1] for v in vic_events]):
#             raise ValueError("Splits events of two demes should be defined relative to one of the demes' id ")
        
#         for i in range(len(vic_events)) :
#             if i == 0 :
#                 continue
#             if vic_events[i][2] in vic_events[i-1][1] and vic_events[i][2] != vic_events[i-1][2]:
#                 raise ValueError("Trying to merge with inactive deme")
            
            
#         # initialize ancestral populations
#         nvic = len(vic_events)
#         ancestrals = [str(vic_events[i][2]) for i in range(nvic)] 
#         count = {}
#         for i, anc in enumerate(ancestrals):
#             cnt = count.get(anc, 0)
#             count[anc] = cnt + 1
#             ancestrals[i] += chr(ord('a') + cnt)

#         ancestrals = ["pop_" + a for a in ancestrals]
#         derived = []
        
        
#         for v in range(nvic):
#             demography_test.add_population(name = ancestrals[v],
#                                            initial_size= com_size[vic_events[v][2]][changetime[vic_events[v][2]].index(vic_dates[v])])
#             tmp = changetime[vic_events[v][2]].index(vic_dates[v]) + 1
#             an_changetime = changetime[vic_events[v][2]][tmp:]
#             an_com_size = com_size[vic_events[v][2]][tmp:]
#             for i in range(len(an_changetime)):
#                 demography_test.add_population_parameters_change(population=ancestrals[v],
#                                                                  time = an_changetime[i],
#                                                                  initial_size= an_com_size[i])
#             derived.extend(vic_events[v][1])

#         #set up split events
#         derived = [str(o) for o in derived]
#         count = {}
#         for i, o in enumerate(derived):
#             cnt = count.get(o, 0)
#             count[o] = cnt + 1
#             if cnt > 0:
#                 derived[i] += chr(ord('a') + cnt - 1) 
#         derived = ["pop_" + d for d in derived]
#         derived = [derived[i*len(derived) // nvic: (i+1)*len(derived) // nvic] for i in range(nvic)]
        
#         for v in range(nvic):
#             demography_test.add_population_split(time = vic_events[v][0], 
#                                                  derived = derived[v], 
#                                                  ancestral = ancestrals[v])
#     demography_test.sort_events()
#     #print(demography_test.debug())


# from collections import Iterable
# def flatten(foo):
#     for x in foo:
#         if hasattr(x, '__iter__') and not isinstance(x, str):
#             for y in flatten(x):
#                 yield y
#         else:
#             yield x

import collections
def flatten(x):
    """
    lizzy need to document this because thx SO
    """
    if isinstance(x, collections.Iterable) and not isinstance(x, str):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

catch = [[1000, [0,1], 1],
          [2000, [1,2], 2]]


flat_list = flatten(catch)
print(flat_list)

print(all([isinstance(x, int) for x in flat_list]))

catch = [["1000", [0,1], 1],
          [2000, [1,2], 2]]

flat_list = flatten(catch)
print(all([isinstance(x, int) for x in flat_list]))

catch = [[1000, [0,1], 1],
          [2000, ["1",2], 2]]

flat_list = flatten(catch)
print(all([isinstance(x, int) for x in flat_list]))

catch = [[1000, [0,1], "1"],
          [2000, [1,2], 2]]

flat_list = flatten(catch)
print(all([isinstance(x, int) for x in flat_list]))