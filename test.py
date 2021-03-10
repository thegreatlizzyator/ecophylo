#!/usr/bin/env python3

from ecophylo import phylogen
from ecophylo import islmodel
from ecophylo import pastdemo
from ecophylo import dosimulate

if __name__ == "__main__":
    import doctest
    print('==================')
    print("phylogen examples")
    print(doctest.testmod(phylogen))
    print('\nphylogen error tests')
    print(doctest.testfile("tests/test-toPhylo.txt"))
    print(doctest.testfile("tests/test-ubranch_mutation.txt"))


    print('\n==================')
    print("islmodel examples")
    print(doctest.testmod(islmodel))
    print('\nislmodel error tests')
    print(doctest.testfile("tests/test-sizes2rates.txt"))
    
    print('\n==================')
    print('Pastdemos examples')
    print(doctest.testmod(pastdemo))
    print('\nPastdemos error tests')
    print(doctest.testfile("tests/test-timeframes.txt"))
    print(doctest.testfile("tests/test-demographic_events.txt"))
    
    print('\n==================')
    print('Dosimulate examples')
    print(doctest.testmod(dosimulate))
    print('\nDosimulate error tests')
    print(doctest.testfile("tests/test-getAbund.txt"))


