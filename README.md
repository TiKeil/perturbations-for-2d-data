# Masterthesis-LOD

## setup

git clone https://github.com/TiKeil/gridlod.git
cd gridlod
git checkout masterthesis
cd ..

git clone https://github.com/TiKeil/Masterthesis-LOD.git

virtualenv -p python2 venv2
. venv2/bin/activate
pip install numpy scipy cython scikit-sparse matplotlib notebook ipython ipdb ipyparallel

if on mac:
brew install suite-sparse
git clone https://github.com/TiKeil/scikit-sparse.git
cd scikit-sparse/
git checkout 0.2-homebrew
python setup.py install
else:
pip install scikits.sparse

echo $PWD > $(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')/gridlod.pth
echo $PWD/Master-thesis-LOD/ > $(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')/Master-thesis-LOD.pth

## work

. venv2/bin/activate
cd Masterthesis-LOD
./start-notebook-server.sh

