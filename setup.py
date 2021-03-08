import setuptools

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

setuptools.setup(
    name="ecophylo",
    version="0.0.6",
    author="Elizabeth Barthelemy",
    author_email="elizabeth.barthelemy@univ-grenoble-alpes.fr",
    url="https://github.com/thegreatlizzyator/ecophylo",
    # download_url="https://github.com/thegreatlizzyator/ecophylo/archive/master.zip",
    description="Coalescent-based simulations of eco-evolutionary biodiverity dynamics",
    long_description=readme,
    long_description_content_type="text/markdown",
    license="CeCILL",
    packages=setuptools.find_packages(),
    install_requires=["setuptools", "requests"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CeCILL-B Free Software License Agreement (CECILL-B)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
