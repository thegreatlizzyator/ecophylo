# EcoPhylo

[![PyPI status](https://img.shields.io/pypi/status/ecophylo.svg)](https://pypi.python.org/pypi/ecophylo/)
[![PyPI version](https://badge.fury.io/py/ecophylo.svg)](https://badge.fury.io/py/ecophylo)
[![PyPI license](https://img.shields.io/pypi/l/ecophylo.svg)](https://pypi.python.org/pypi/ecophylo/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/ecophylo.svg)](https://pypi.python.org/pypi/ecophylo/)

project description

The source for this project is available [here][src].

Provide a link to full documentation tutorial [here][tutorial]. 

----

Publications using this package here. 

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

### Installation and usage in R (on Windows aka the painfull way)

> For François : First create a new project in rstudio with **version control** from the [github](https://github.com/thegreatlizzyator/ecophylo). 

You will need to install `{reticulate}` package and python dependencies.

```r
# NOT RUN
install.packages('reticulate') # answer yes
library(reticulate)
# python dependencies
conda_install('r-reticulate', c('msprime','ete3','pandas'))
```

## Running ##

### In python

In python (version >= 3.7) you can simulate trees using this package main function in python. 
```python
from ecophylo import dosimuls
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True)
```

Note that there is random events during the simulation that you can control by setting a seed :
```python
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True, seed = 42)
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True, seed = 42)
```

### In R

**A vignette is being written at this moment.**
