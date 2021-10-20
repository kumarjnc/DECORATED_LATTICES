# Summary 
## Code to calculate dispersion, density of state (DOS), local density of state (LDOS) of decorated square lattices solving tight binding models.

### Installation
```bash
git clone https://github.com/kumarjnc/DECORATED_LATTICES.git
```

# Example
## Band dispersion of Decorated lattice
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
# Visualising dispersion, DOS and LDOS
# Reshaping the output
```python
dispersion_data=np.loadtxt("dispersion_lieb.dat")
kvecs = np.loadtxt("kvecs.dat")
X,Y=np.meshgrid(kvecs[:,0],kvecs[:,1])
Zup=dispersion_data[:,0].reshape(Nk,Nk)
Zdw=dispersion_data[:,1].reshape(Nk,Nk)
Zflat=dispersion_data[:,2].reshape(Nk,Nk)
Ztop=dispersion_data[:,3].reshape(Nk,Nk)
```
# 3d Plot interactive
```python
%matplotlib notebook
from matplotlib import interactive
interactive(True)
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.
ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('$E$')


surf = ax.plot_surface(X, Y,Zup, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
surf = ax.plot_surface(X, Y,Zdw, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
surf = ax.plot_surface(X, Y,Zflat, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
surf = ax.plot_surface(X, Y,Ztop, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()
```

