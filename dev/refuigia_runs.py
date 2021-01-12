# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:54:00 2020

@author: barthele
"""
import all_ecophylo

times = [50, 100, 200, 500, 1200, 1500, 2000, 3000, 5000, 10000]

for t in times:
    names = f'refugia_trees_mu_{t}.txt'
    print(names)
    all_ecophylo.dosimuls(1000, 500, [10000], [0.001], pastprior= [1000, 5000], changetime = [t], savetrees = True, saveto = names)
    

test = all_ecophylo.dosimuls(100, 100, [10000], [0.001], savetrees = False, saveto = "")