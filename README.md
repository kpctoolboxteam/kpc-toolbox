# KPC-Toolbox
KPC-Toolbox: MATLAB and Python toolbox to fit Markovian Arrival Processes

Current version: 0.5.0

The toolbox is now dual-language: the reference MATLAB implementation lives in
`matlab/`, and a standalone Python package lives in `python/` (see
`python/README.md`).

The Python package is an AI-assisted port of the MATLAB toolbox; see
AI-DISCLAIMER.md.

Website: http://www.cs.wm.edu/MAPQN/kpctoolbox.html

This software is released under the BSD-3 license, see LICENSE.TXT.

If you are using the KPC-Toolbox for scientific papers or technical reports, please consider citing the following publications:

[1] G.Casale, E.Z.Zhang, E.Smirni. 
KPC-Toolbox: Best Recipes for Automatic Trace Fitting Using Markovian Arrival Processes 
Elsevier Performance Evaluation, 67(9):873-896, Sep 2010.

[2] G.Casale, E.Z.Zhang, E.Smirni. 
Trace Data Characterization and Fitting for Markov Modeling
Elsevier Performance Evaluation, 67(2):61-79, Feb 2010.

GETTING STARTED (MATLAB)

Add all directories to the MATLAB classpath, for example using the command:

```
addpath(genpath('MY_INSTALLATION_PATH/kpc-toolbox/matlab'))
```

A demonstrator of the tool can be run using the command

```
demo_run
```

To get help for the Markovian arrival process fitting tool, type

```
help kpcfit_auto
```

For the phase-type distribution fitting tool, type

```
help kpcfit_ph_auto
```

GETTING STARTED (PYTHON)

Install the package from PyPI, which also brings in line-solver, numpy and scipy:

```
pip install kpc-toolbox
```

A demonstrator of the tool can be run using the commands

```python
from kpctoolbox import demo
demo.demo_run()
```

To fit a Markovian arrival process to a trace, use

```python
import kpctoolbox as k
from kpctoolbox import demo

S = demo.load_trace('BCAUG89')               # interarrival-time samples
trace = k.kpcfit_init(S)
MAP = k.kpcfit_auto(trace, OnlyAC=True)[0]   # best fit, as (D0, D1)
```

and for the phase-type distribution fitting tool, use

```python
import numpy as np
import kpctoolbox as k

D0, D1 = k.map_hyperexp(1.0, 8.0, 0.6)                            # target
E = np.array([float(k.map_moment(D0, D1, i)) for i in (1, 2, 3)])
fits = k.kpcfit_ph_auto(E, k.kpcfit_ph_options(E))                # ranked fits
PH = fits[0][0]                                                   # best (D0, D1)
```

The Python API notes, including the few deliberate differences from the MATLAB
API, are in `python/README.md`.
