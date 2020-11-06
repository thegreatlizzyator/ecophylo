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
or install from github :
```shell
python3 -m pip install git+https://github.com/thegreatlizzyator/ecophylo/tree/packaging/dist/ecophylo-0.0.5.tar.gz
```

## Running ##

In python ( version >= 3.6) you can simulate trees using this package main function :

```python
from ecophylo import dosimuls
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True)
```

Note that there is random event during the simulation that you can control by setting a seed :
```python
np.random.seed(42)
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True)
np.random.seed(42)
dosimuls(nsim = 5, sample_size = 100, comprior = [1000,10e9], muprior = [1e-6] , verbose = True)
```
