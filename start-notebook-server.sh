#!/bin/bash
#
# ~~~
# This file is part of the project: perturbations-for-2d-data
#   https://github.com/TiKeil/perturbations-for-2d-data.git
# Copyright 2017-2020: Tim Keil. All rights reserved.
# License: Licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Tim Keil (2017, 2020)
# ~~~

jupyter-notebook --ip 0.0.0.0 --notebook-dir=$PWD/notebooks/ --port=18081
