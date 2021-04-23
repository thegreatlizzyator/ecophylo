#!/usr/bin/env python3.7

from ecophylo import phylogen
from ecophylo import pastdemo
from ecophylo import dosimulate

if __name__ == "__main__":
    import doctest
    print('==================')
    print("phylogen examples")
    print(doctest.testmod(phylogen))
    print('\nphylogen error tests')
    print('toPhylo :', doctest.testfile("tests/test-toPhylo.txt"))
    print('ubranch_mutation :', doctest.testfile("tests/test-ubranch_mutation.txt"))


    print('\n==================')
    print('Pastdemos examples')
    print(doctest.testmod(pastdemo))
    print('\nPastdemos error tests')
    print('timeframes :', doctest.testfile("tests/test-timeframes.txt"))


    print('\n==================')
    print('Dosimulate examples')
    print(doctest.testmod(dosimulate))
    print('\nDosimulate error tests')
    print('simulate :', doctest.testfile("tests/test-simulate.txt"))
    print('getAbund :', doctest.testfile("tests/test-getAbund.txt"))


