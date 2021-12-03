import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset

with Dataset("nc/isla.nc","r") as ifile:
    η0=ifile["η0"][:]
    h0=ifile["h0"][:]
    xi=ifile["x"][:]

#cantidad de puntos
kcells=len(h0)
h=h0+η0
dx=5
dt=0.01
u=np.zeros((3,kcells+1))
eta=np.zeros((3,kcells))
eta[0]=η0
g=9.81
ku=2*g*dt/dx
kh=2*dt/dx
eps=0.03
hmin=0.15
mu=dt*np.sqrt(g*h0)/dx
print('mu:', mu)

plt.figure(figsize=(9,4))
#pasos de tiempo
for n in np.arange(500.5/dt):
    ni=int(n)%3
    if n==0:
        #celdas parcialmente secas
        idx_dry=np.asarray(h<hmin).nonzero()[0]
        print(n, 'idx', idx_dry)
        lim=(i for i in idx_dry if (h[i-1]>=hmin or h[i+1]>=hmin) )
        lim = list(lim)
        dry = ( i for i in idx_dry if i not in lim)
        dry = list(dry)
        print('lim:', lim, 'dry:', dry)
    else:
    #if n!=0:
        if n==1:
            u[1][1:-1]=u[0][1:-1]-ku*np.diff(eta[0])/2
            hu=h*u[0][1:]
            hu=np.insert(hu,0, 0)
            eta[1]=eta[0]-kh*np.diff(hu)/2
        elif n>1:
            u[ni][1:-1]=u[ni-2][1:-1]-ku*np.diff(eta[ni-2])
            hu=h*u[ni-2][1:]
            hu=np.insert(hu,0, 0)
            eta[ni]=eta[ni-2]-kh*np.diff(hu)
        u[ni][0]=0#u[ni][1]
        u[ni][-1]=0#u[ni][-2]
        h=eta[ni]+h0
        #celdas parcialmente secas
        idx_dry=np.asarray(h<hmin).nonzero()[0]
        print(n, 'idx', idx_dry)
        lim=(i for i in idx_dry if (h[i-1]>=hmin or h[i+1]>=hmin) )
        lim = list(lim)
        dry = ( i for i in idx_dry if i not in lim)
        dry = list(dry)
        print('lim:', lim, 'dry:', dry)
        ucond=u[ni][1:]
        ucond[h<=hmin]=0
        u[ni][1:]=ucond
            
        #filtrado
        eta[ni][1:-1]=(1-eps)*eta[ni][1:-1]+0.5*eps*(np.roll(eta[ni],1)[1:-1]+np.roll(eta[ni],-1)[1:-1])
        eta[ni][dry]=-h0[dry]
        h=eta[ni]+h0

    #graficación
    if n%100==0:
        #crea eta mas adecuada para graficación
        eta_plot=np.copy(eta[ni])
        eta_plot[dry]=np.NaN
        plt.plot(xi, -h0, 'b')
        plt.plot(xi, eta_plot,'.r')
        plt.axis([0, kcells*dx, -1, 1.5])
        plt.title('t={}s, eps={}'.format(n*dt, eps))
        plt.plot(xi, u[ni][1:],'k')
        if n==0:
            plt.pause(1.5)
        else:
            plt.pause(0.1)
        plt.cla()
plt.show()

