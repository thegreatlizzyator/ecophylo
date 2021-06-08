rm dist/*
python3 setup.py sdist bdist_wheel
printf "\n=========================================================================\n\n"
python3 -m twine check dist/*
printf "\n=========================================================================\n"
printf "=============================  Install   ================================\n"
printf "=========================================================================\n"
printf "installing without internet. will crash if change in dependencies\n"
python3.7 -m pip install dist/ecophylo-0.1.3.tar.gz --no-deps --no-index
#python3 -m pip install dist/ecophylo-0.1.3.tar.gz 
printf "\n=========================================================================\n"
printf "=============================  TESTING   ================================\n"
printf "=========================================================================\n"
python3.7 test.py

# python3.7 -m sphinx.cmd.build -b html /home/maxime/BEE/ecophylo/docs /home/maxime/BEE/ecophylo/docs/_build
# from https://brendanhasz.github.io/2019/01/05/sphinx.html#cross-referencing
