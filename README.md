# EcoPhylo

project description

The source for this project is available [here][src].

Provide a link to full documentation tutorial [here][tutorial]. 

----

Publications using this package here. 

[src]: https://github.com/thegreatlizzyator/ecophylo

## Installation ##

To install this package, you can either download the tar.gz file in the dist directory of the repository :
```shell
python3 -m pip install <path_to>/ecophylo-<VERSION>.tar.gz
``` 

<!---
TODO : put the package on Pypi
or instal from the Python Package index (Pypi) :

```shell
python3 -m pip install ecophylo
``` 
-->

or install from github (it does not work for private repository):
```shell
python3 -m pip install git+https://github.com/thegreatlizzyator/ecophylo/tree/packaging/dist/ecophylo-0.0.5.tar.gz
```
<!---
TODO : check if better way to check on dependencies
-->
To end the installation please check for the dependancies. In python you can check the dependencies like this :
```python
# dependencies
import msprime
import numpy as np
import sys
from ete3 import Tree
import pandas as pd
```

You can install them by running (numpy is an example):
```shell
python3 -m pip install numpy
```
You can experiment difficulties for installing msprime, check the help [here][msprime]

[msprime]: https://msprime.readthedocs.io/en/stable/installation.html


## Running ##

In python ( version >= 3.6) you can simulate trees using this package main function 
```python
from ecophylo import dosimuls
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True)
```

Note that there is random events during the simulation that you can control by setting a seed :
```python
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True, seed = 42)
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True, seed = 42)
```
