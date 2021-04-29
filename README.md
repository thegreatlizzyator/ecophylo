# EcoPhylo

[![PyPI status](https://img.shields.io/pypi/status/ecophylo.svg)](https://pypi.python.org/pypi/ecophylo/)
[![PyPI version](https://badge.fury.io/py/ecophylo.svg)](https://badge.fury.io/py/ecophylo)
[![PyPI license](https://img.shields.io/pypi/l/ecophylo.svg)](https://pypi.python.org/pypi/ecophylo/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/ecophylo.svg)](https://pypi.python.org/pypi/ecophylo/)

*Ecophylo* is dedicated to coalescent-based simulation of eco-evolutionary dynamics. Species assemblages and their shared ancestry can be simulated by jointly taking into account the influence of past demographic fluctuations and extinctions along with how divergent genotypes have introduced new species over time through speciation.

The source for this project is available [here][src].

----

[src]: https://github.com/thegreatlizzyator/ecophylo

## Dependencies

This package depend on **python 3.7**.

There are multiple dependencies to other python packages. This will be done when installing ecophylo tool in most cases. 

You can experiment difficulties for installing msprime, check the help [here][msprime]

[msprime]: https://tskit.dev/msprime/docs/stable/installation.html

## Installation (with pip)

To install this package, you can either download the tar.gz file in the *dist* directory of the repository :
```shell
python3 -m pip install <path_to>/ecophylo-<VERSION>.tar.gz
``` 

or install from the Python Package index (Pypi) :

```shell
python3 -m pip install ecophylo
``` 

or install from github (it does not work for private repository):

```shell
python3 -m pip install git+https://github.com/thegreatlizzyator/ecophylo/tree/packaging/dist/ecophylo-0.0.5.tar.gz
```

## Running ##

You can simulate trees using this package main function in python.

```python
import ecophylo
n = 25 #the number of sampled individuals
com_size = [[5000,10000,50000]] # the size of the assemblage in the past, the first element is the current assemblage size
mu = 0.001 # the point mutation rate
changetime = [[0,700,10000]] # the dates (in generation time) at which the assemblage has changed sizes in the past


tree = ecophylo.simulate(samples = n,
                         com_size= com_size, 
                         mu = mu, 
                         changetime= changetime, seed= 42)

print(tree)
```
