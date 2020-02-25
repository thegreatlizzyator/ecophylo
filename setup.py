# -*- coding: utf-8 -*-
"""
find cool name.

"The unpredictable and the predetermined unfold together to make everything the way it is.
It’s how nature creates itself, on every scale, the snowflake and the snowstorm."
Tom Stoppard, Arcadia (Faber and Faber, London, 1993)

@author: Elizabeth Barthelemy
"""

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

setup(name='findcoolname',
      version='0.1',
      description='Description',
      url='http://github.com/thegreatlizzyator/findcoolname',
      author='Elizabeth Barthelemy',
      author_email='elizabeth.barthelemy@univ-grenoble-alpes.fr',
      license='CeCILL',
      packages=find_packages(exclude=('tests', 'docs')),
      zip_safe=False)