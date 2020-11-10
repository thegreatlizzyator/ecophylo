#!/usr/bin/env python3

from ecophylo import islmodel
from ecophylo import pastdemo

if __name__ == "__main__":
    import doctest
    print("Testing islmodel functions")
    print(doctest.testmod(islmodel))
    print("Testing pastdem functions")
    print(doctest.testmod(pastdemo))

