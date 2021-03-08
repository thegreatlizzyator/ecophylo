rm dist/*
python3 setup.py sdist bdist_wheel
python3 -m pip install dist/ecophylo-0.0.6.tar.gz
python3 test.py
