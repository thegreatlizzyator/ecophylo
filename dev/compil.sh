rm dist/*
python3 setup.py sdist bdist_wheel
python3 -m pip install dist/ecophylo-0.0.7.tar.gz
printf "\n=========================================================================\n"
printf "=============================  TESTING   ================================\n"
printf "=========================================================================\n"
python3 test.py
