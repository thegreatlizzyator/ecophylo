from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

setup(
    name="ecophylo",
    version="0.0.1",
    url="https://github.com/thegreatlizzyator/ecophylo",
    download_url="https://github.com/thegreatlizzyator/ecophylo/archive/master.zip",
    license="CeCILL",
    author="Elizabeth Barthelemy",
    author_email="elizabeth.barthelemy@univ-grenoble-alpes.fr",
    description="Coalescent-based simulations of eco-evolutionary biodiverity dynamics",
    long_description=readme,,
    packages=find_packages(exclude=("tests",)),
    install_requires=["setuptools", "requests"],
    classifiers=[
	"Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: CeCILL FREE SOFTWARE LICENSE AGREEMENT",
    ],
)