# Summary 
## Code to calculate dispersion, density of state (DOS), local density of state (LDOS) of decorated lattices solving tight binding models

### Installation
```bash
git clone https://github.com/kumarjnc/DECORATED_LATTICES.git
```

# Example
## Band dispersion of Lieb lattice
```python
import numpy as np 
from numpy import tile
import matplotlib.pyplot as plt
from scipy import linalg as LA
import math,cmath
from scipy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
exp=np.exp
sqrt = math.sqrt
log=np.log
from numpy.linalg import inv
from numpy import zeros
import main
import params
main.band(params.t1,params.t2,params.a,params.Nk,params.density,params.mass,params.dim,params.totomega,params.m,params.sigma,square=True)
``` 
