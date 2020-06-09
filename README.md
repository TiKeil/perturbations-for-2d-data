[![DOI](https://zenodo.org/badge/270818299.svg)](https://zenodo.org/badge/latestdoi/270818299)

```
# ~~~
# This file is part of the project: perturbations-for-2d-data
#   https://github.com/TiKeil/perturbations-for-2d-data.git
# Copyright 2017-2020: Tim Keil. All rights reserved.
# License: Licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Tim Keil (2017 - 2020)
#
# ~~~
```

# perturbations-for-2d-data
This repository contains python code for perturbing a two-dimensional image. The code originally belongs to the master thesis "Variational crimes in the Localized orthogonal decomposition method" by Tim Keil. 
We provide some small tutorials and examples 
The entire idea for two dimensional data is based on the python module 'gridlod'. This module has been used and developed by Fredrik Hellman on https://github.com/fredrikhellman/gridlod. Our repository provides classes for building data for gridlod but can be used independently for arbitrary two dimensional data and images. 

## Setup

We recommend to use a virtual python3 environment, clone our repository and simply run

```
pip install -e . 
pip install -r requirements.txt
```

and start a notebook servia via

```
./start-notebook-server.sh
```
