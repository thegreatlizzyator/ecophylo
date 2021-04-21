#!/usr/bin/env python3.7

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
    print('toPhylo :', doctest.testfile("tests/test-toPhylo.txt"))
    print('ubranch_mutation :', doctest.testfile("tests/test-ubranch_mutation.txt"))


    # print('\n==================')
    # print("islmodel examples")
    # print(doctest.testmod(islmodel))
    # print('\nislmodel error tests')
    # print('sizes2rates :', doctest.testfile("tests/test-sizes2rates.txt"))
    # print('mergesizes2rates :', doctest.testfile("tests/test-mergesizes2rates.txt"))
    # print('population_configurations_stripe :', doctest.testfile("tests/test-population_configurations_stripe.txt"))


    print('\n==================')
    print('Pastdemos examples')
    print(doctest.testmod(pastdemo))
    print('\nPastdemos error tests')
    print('timeframes :', doctest.testfile("tests/test-timeframes.txt"))
    print('demographic_events :', doctest.testfile("tests/test-demographic_events.txt"))


    print('\n==================')
    print('Dosimulate examples')
    print(doctest.testmod(dosimulate))
    print('\nDosimulate error tests')
    print('simulate :', doctest.testfile("tests/test-simulate.txt"))
    print('getAbund :', doctest.testfile("tests/test-getAbund.txt"))


