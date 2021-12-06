'''
Ejemplo de cómo crear un netCDF a partir de un arreglo

crea batimetría para ejemplo de ladera
'''

import numpy as np
from netCDF4 import Dataset

#define arreglos a almacenar

#cantidad de puntos
kcells=201
#define nivel estable
lvl0=10
#define pendiente
step=10/201
#define batimetría
bat=np.array(list(np.arange(20,10,-step)))
bat[90:111]=14
bat=-bat
print('bat', bat.shape)
#define η0
η0=np.zeros(bat.shape)
η0[0:40]=1.0
η0-=np.minimum(0,bat)

dx=5
x= kcells*dx
xi=np.arange(0,x,dx)
print(η0.shape, bat.shape, xi.shape)
with Dataset("ladera.nc", "w") as root:
    root.description="batimetría para ejemplo de isla"
    dim_k= root.createDimension("kcells",kcells)
    var_η0= root.createVariable("η0","f8",("kcells",))
    var_h0= root.createVariable("h0","f8",("kcells",))
    var_x= root.createVariable("x","f8",("kcells",))
    var_η0[:]=η0
    var_h0[:]=bat
    var_x[:]=xi
