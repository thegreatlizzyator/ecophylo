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

# def dosimuls(nsim, sample_size, comprior, muprior, lim_mrca = None, sstype="SFS",
#              prior_distrib = "uniform", npop=1, withmigr=False, init_ratesprior=None,
#              init_sizeprior=None, pastprior=None, changetime = None,
#              nsplit=None, massprior=None, migrfrom=None, migrto=None,
#              verbose=False, savetrees= False, saveto = "", seed = None):
#     pass


import ecophylo as eco

def dosimuls(nsim, samples, com_size, mu, init_rates = None, changetime = None,
             mrca = None, migr = None, migr_time = None,
             vic_events = None, 
             verbose = False, sumstat = None, result = ['Params'], 
             file = None, seed = None):
    


    samples, com_size, mu, init_rates, changetime, mrca, migr, migr_time, \
                    vic_events, verbose, seed, prior_locate = eco.check_params(
        samples = samples, com_size = com_size, mu = mu, init_rates = init_rates, 
        changetime = changetime, mrca = mrca, migr = migr, 
        migr_time = migr_time, vic_events = vic_events, 
        verbose = verbose, seed = seed, prior_locate = 'naive'
    )
    i = 1
    print(com_size, prior_locate)
    l = []
    l.append(com_size[prior_locate[0][1]][prior_locate[0][2]])
    while i <= nsim:
        samples, com_size, mu, init_rates, changetime, mrca, migr, migr_time, \
                    vic_events, verbose, seed, prior_locate = eco.check_params(
            samples = samples, com_size = com_size, mu = mu, init_rates = init_rates, 
            changetime = changetime, mrca = mrca, migr = migr, 
            migr_time = migr_time, vic_events = vic_events, 
            verbose = verbose, seed = seed, prior_locate = prior_locate
        )   
        l.append(com_size[prior_locate[0][1]][prior_locate[0][2]])
        i+=1
    # l = "Fini"
    return l

print(dosimuls(nsim = 5, samples = [10, 10],
    com_size = [[500, 1000], [200, [150, 500, "uniform"], 600]],
    mu = 0.001,
    changetime = [[0, 200], [0, 100, 500]],
    seed = 42))

print(dosimuls(nsim = 5, samples = [10, 10],
    com_size = [[500, 1000], [200, [150, 500, "uniform"], 600]],
    mu = 0.001,
    changetime = [[0, 200], [0, 100, 500]]))

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

print()
a = 2
def foo(a):
    a = a +1
    print(a)

foo(a = a)
print(a)

