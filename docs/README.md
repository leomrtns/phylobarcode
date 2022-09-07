<img src="phylobarcode.png" height="150" align="right">

# phylobarcode

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Introduction
**Phylobarcode** is a tool to search and analyse long operons with phylogenetic signal.

It can search for potential long segments ('operons' or 'barcodes') which can serve as phylogenetic markers.
It can then search for potential primers which can be used to amplify the operon.
And it can be used to reconstruct the evolutionary history of the amplicons.


## Installation


### Requirements

* `conda`
* linux
* python=3.9

The other requirements can be installed with `conda`. In particular `parasail-python` is not available yet for python 3.10

### Generate a conda environment

This software depends on several other packages, installable through conda or pip.
The suggested installation procedure is to create a conda environment (to take care of dependencies) and then installing
the python package:
```bash
conda update -n base -c defaults conda # probably not needed, but some machines complained about it
conda env create -f environment.yml  
conda activate phylobarcode
python setup.py install # or "pip install ." 
```

Since this software is still under development, these two commands are quite useful:
```bash
conda env update -f environment.yml # update conda evironment after changing dependencies
# installs in development mode (modifications to python files are live):
pip install -e .  # or python setup.py develop
```

## License 
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2022-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

phylobarcode is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).
