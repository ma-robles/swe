'''
Ejemplo de cómo crear un netCDF a partir de un arreglo

crea batimetría para ejemplo de isla 1D
'''

import numpy as np
from netCDF4 import Dataset

#define arreglos a almacenar

#cantidad de puntos
kcells=201
#define nivel estable
lvl0=10
#define batimetría
bat =np.ones(kcells)*lvl0
#define pendiente
step=10.5/101
for i,g in enumerate(bat):
    if i<=100:
        bat[i+1]=g-step
    elif i<200:
        bat[i+1]=g+step
#define η0
η0=np.zeros(bat.shape)
η0[0:40]=1.0
η0-=np.minimum(0,bat)

dx=5
x= kcells*dx
xi=np.arange(0,x,dx)
print(η0.shape, bat.shape, xi.shape)
with Dataset("isla.nc", "w") as root:
    root.description="batimetría para ejemplo de isla"
    dim_k= root.createDimension("kcells",kcells)
    var_η0= root.createVariable("η0","f8",("kcells",))
    var_h0= root.createVariable("h0","f8",("kcells",))
    var_x= root.createVariable("x","f8",("kcells",))
    var_η0[:]=η0
    var_h0[:]=bat
    var_x[:]=xi
