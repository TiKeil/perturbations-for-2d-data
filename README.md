# Masterthesis-LOD

```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
```

This repository contains the python coding that belongs to the master thesis "Variational crimes in the Localized orthogonal decomposition method" by Tim Keil. Furthermore, it presents the generation of each python figure that is included in the thesis. If you intend to try the tests and the generation of the figures by yourself you need to follow the setup steps that are explained in the next section. After the setup, you can either run the jupyter-notebooks or run the complete python scripts by yourself. Both strategies allow for adjusting the test to own purposes such as different coefficients or perturbations.
The entire coding is based on the python module 'gridlod'. This module has been used and developed by Axel Målqvist and Fredrik Hellmann. Our repository provides further classes, that work as an extension of 'gridlod' and enable an efficient construction of two dimensional coefficients. 

## Setup

First you need to select your favorite folder and clone the git repository that contains the 'gridlod' module. Furthermore, you switch to the compatible branch on which the extensions are established.

```
git clone https://github.com/TiKeil/gridlod.git
cd gridlod
git checkout masterthesis
cd ..
```

Now, add our repository with
``` 
git clone https://github.com/TiKeil/Masterthesis-LOD.git
```

In order to connect 'gridlod' to our new files, you need to work with a virtual environment. First, you construct this environment and then you activate it.

```
virtualenv -p python2 venv2
. venv2/bin/activate
```

In this 'virtualenv', you have to install the required packages for python and everything that remains to run 'gridlod' and our files. 

```
pip install numpy scipy cython scikit-sparse matplotlib notebook ipython ipdb ipyparallel
```

If you are working on a mac, you need to use a slightly hacked version of scikit-sparse and do

```
brew install suite-sparse
git clone https://github.com/TiKeil/scikit-sparse.git
cd scikit-sparse/
git checkout 0.2-homebrew
python setup.py install
```

On Linux, just install scikits.sparse

```
pip install scikits.sparse
```

Now, link gridlod to our folder that enables to use it as a module.

```
echo $PWD > $(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')/gridlod.pth
echo $PWD/Master-thesis-LOD/ > $(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')/Master-thesis-LOD.pth
```

## Work

Every time you want to work with the virtualenv, you simply need to run 

```
. venv2/bin/activate
cd Masterthesis-LOD
```

We recommend to work with the jupyter-notebooks. They provide a good overview of the content and simplify manual changes.  

```
./start-notebook-server.sh
```

Just copy the url into your favorite Browser, click on 'notebooks' and have a look at the prepared files. The interface enables a better insight to the, coding due to further explanations.

Clearly, you can also run the complete python scripts in the folder 'virtualenv_files' if you do not want to use jupyter-notebooks. For example, you can do the following commands with activated 'virtualenv'

```
cd virtualenv_files
cd generate_figures
python 2.4.1_MsExampleFEM1d
```
