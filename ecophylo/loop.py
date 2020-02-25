# -*- coding: utf-8 -*-
#
# This file is a part of the ecophylo package.
# 
# ecophylo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the CeCILL FREE SOFTWARE LICENSE AGREEMENT
# along with ecophylo.

"""
Module responsible for simulating phylogenies along braod ranges of parameter
values that can be used later for retreiving the distribution of most likely
parameter values using Approximate Bayesian Computation.

@author: Elizabeth Barthelemy
 """
import numpy as np
import pandas as pd
import sys
import os
import datetime
import h5py

import simulate
import sumstats


def dosimuls(nsim, sample_size, comprior, muprior, sstype="SFS", npop=1,
             nepoch=1, withmigr=False, init_ratesprior=None,
             init_sizeprior=None, pastprior=None, maxtime=None,
             nsplit=None, massprior=None, migrfrom=None, migrto=None,
             verbose=False, writehdf5 = False, filename ='', metadata =''):

    # CHECKS HERE FOR IDIOT-PROOFING

    # GENERATING PRIOR DISTRIBUTION TABLE
    comments = 'Simulating eco-evolutionary dynamics over the following parameters ranges:'
    params.count = 0
    df = pd.DataFrame()

    if muprior is not None:
        df['mu'] = params(muprior, nsim)
        comments += (f'\n -mutation rate in {muprior}')

    if comprior is not None:
        df['comsize'] = params(comprior, nsim)
        comments += (f'\n -community size in {comprior}')

    if pastprior is not None:
        comments += (f'\n -{nepoch} epoch(s) in which community size varies in {pastprior}')
        pastsizes = [params(pastprior, nsim) for _ in range(nepoch)]
        colnames = ['pastsize{}'.format(i) for i in range(1, nepoch+1)]
        for i in range(nepoch):
            df[colnames[i]] = pastsizes[i]

    if init_sizeprior is not None and init_ratesprior is not None and npop is not None:
        initrates = [params(init_ratesprior, nsim) for _ in range(npop)]
        colnames = ['initrate{}'.format(i) for i in range(1, npop+1)]
        for i in range(npop):
            df[colnames[i]] = initrates[i]

        initsizes = [params(init_sizeprior, nsim) for _ in range(npop)]
        colnames = ['initsize{}'.format(i) for i in range(1, npop+1)]
        for i in range(npop):
            df[colnames[i]] = initsizes[i]
        comments += (f'\n -{npop} population(s) in which initial community sizes vary in {init_sizeprior} and growth rates vary in {init_ratesprior}')

    if withmigr and npop is not None:
        df['m'] = params([0, 1], nsim)
        comments += ('\n -there is migration between populations')

    if massprior is not None:
        massdates = [params(massprior, nsim) for _ in range(nsplit)]
        colnames = ['splitdate{}'.format(i) for i in range(1, nsplit+1)]
        for i in range(nsplit):
            df[colnames[i]] = massdates[i]
        comments += (f'\n -populations {migrfrom} and {migrto} split at dates ranging in {massprior}')

    comments += f"\nEstimating the associated model's {params.count} parameter(s) based on the {sstype} summary statistic(s)"
    if verbose:
        print(comments)

    # LOOP OVER PRIOR DISTRIBUTIONS AND SIMULATE PHYLOGENIES
    ss = list()
    for i in range(nsim):
        mu = df.loc[i, 'mu']
        com_size = df.loc[i, 'comsize']
        try:
            m = df.loc[i, 'm']
        except KeyError:
            m = None
        past_sizes = list(df[[col for col in df if col.startswith('pastsize')]].iloc[i, ])
        init_rates = list(df[[col for col in df if col.startswith('initrate')]].iloc[i, ])
        init_sizes = list(df[[col for col in df if col.startswith('initsize')]].iloc[i, ])
        split_dates = list(df[[col for col in df if col.startswith('splitdate')]].iloc[i, ])

        phylo = simulate.simulate(sample_size, com_size, mu, npop, nepoch, m, init_rates, init_sizes, past_sizes, maxtime, split_dates, migrfrom, migrto, seed=1, verbose=False)

        # compute summary statistics (only SFS is supported for now)
        if sstype == 'SFS':
            ss.append(sumstats.getSFS(phylo))
    if sstype == 'SFS':
        SFS = np.zeros((nsim, max(len(x) for x in ss)))
        for i, j in enumerate(ss):
            SFS[i][0:len(j)] = j
    ssdf = pd.DataFrame(SFS)

    # export to hdf5 format if specified
    if writehdf5:
        save2path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Simulations")
        if not os.path.exists(save2path):
            os.makedirs(save2path, exist_ok=True)
        filename = os.path.join(save2path,f'{filename}ecophylo_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S" +'.h5'))
        hdf_file = h5py.File(filename, 'w')
        hdf_file.create_dataset('params', data = df)
        hdf_file.create_dataset('sumstats', data = ssdf)
        hdf_file.close()
        
    return df, ssdf

def counter(func):
    func.count += 1

def params(lim, nsim, distrib = "uniform"):
    """
    Make a list of parameters to test from prior distribution limits.

    Parameters
    ----------
    lim: list
        the inferior and superior limits of the prior distribution.
    nsim: int
        number of replicates.
    distrib: str, optional
        the type of prior distribution. The default is "uniform".

    Returns
    -------
    params: list
        a list of nsim parameter values.

    """
    if len(lim) != 2:
        sys.exit('inferior and superior prior limits are needed')
    counter(params)
    #only supports uniform prior distributions at the moment
    if distrib == "uniform":
        p = np.random.uniform(lim[0], lim[1], size = nsim)
    return p