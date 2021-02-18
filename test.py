#!/usr/bin/env python3

from ecophylo import toPhylo
from ecophylo import islmodel
from ecophylo import pastdemo
from ecophylo import dosimulate

if __name__ == "__main__":
    import doctest
    print("Testing toPhylo functions")
    print(doctest.testmod(toPhylo))
    print("\n\nTesting islmodel functions")
    print(doctest.testmod(islmodel))
    print("\n\nTesting pastdem functions")
    print(doctest.testmod(pastdemo))
    print("\n\nTesting dosimulate functions")
    print(doctest.testmod(dosimulate))

