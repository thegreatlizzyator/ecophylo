import setuptools

# import pathlib

# here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
# long_description = (here / 'README.md').read_text(encoding='utf-8')

# with open("README.md", "r") as readme_file:
#     readme = readme_file.read()
import os
current_path = os.path.abspath(os.path.dirname(__file__))
def read_file(*parts):
    with open(os.path.join(current_path, *parts), encoding='utf-8') as reader:
        return reader.read()

setuptools.setup(
    name="ecophylo",
    version="0.1.2",
    author="Elizabeth Barthelemy",
    author_email="elizabeth.barthelemy@univ-grenoble-alpes.fr",
    url="https://github.com/thegreatlizzyator/ecophylo",
    # download_url="https://github.com/thegreatlizzyator/ecophylo/archive/master.zip",
    description="Coalescent-based simulations of eco-evolutionary biodiverity dynamics",
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    license="CeCILL",
    packages=setuptools.find_packages(),
    install_requires=["setuptools", "requests"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.7",
        "License :: CeCILL-B Free Software License Agreement (CECILL-B)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
