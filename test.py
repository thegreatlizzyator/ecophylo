#!/usr/bin/env python3.7

from ecophylo import phylogen
from ecophylo import pastdemo
from ecophylo import dosimulate
from ecophylo import sumstat

if __name__ == "__main__":
    import doctest
    print('==================')
    print("Phylogen examples")
    print(doctest.testmod(phylogen))
    print('* Phylogen error tests')
    print('toPhylo :', doctest.testfile("tests/test-toPhylo.txt"))
    print('ubranch_mutation :', doctest.testfile("tests/test-ubranch_mutation.txt"))


    print('\n==================')
    print('Pastdemos examples')
    print(doctest.testmod(pastdemo))
    print('* Pastdemos error tests')
    print('timeframes :', doctest.testfile("tests/test-timeframes.txt"))

    print('\n==================')
    print('Dosimulate examples')
    print(doctest.testmod(dosimulate))
    print('* Dosimulate error tests')
    print('simulate :', doctest.testfile("tests/test-simulate.txt"))
    print('simulate :', doctest.testfile("tests/test-dosimuls.txt"))
    print('simulate :', doctest.testfile("tests/test-sample.txt"))

    print('\n==================')
    print('Sumstat examples')
    print(doctest.testmod(sumstat))
    print('* Sumstat error tests')
    print('getAbund :', doctest.testfile("tests/test-getAbund.txt"))
    print('getDeme :', doctest.testfile("tests/test-getDeme.txt"))

