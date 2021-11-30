import numpy as np
from matplotlib import pyplot as plt

#cantidad de puntos
kcells=101
#define nivel estable
lvl0=10
#define tierra
bat =np.ones(kcells)*lvl0
step=10.5/51
for i,g in enumerate(bat):
    if i<=50:
        bat[i+1]=g-step
    elif i<100:
        bat[i+1]=g+step
#define h0
h0=np.copy(bat)
#h0=np.ones(bat.shape)*10
eta0=np.zeros(h0.shape)
eta0[0:20]=1.0
eta0-=np.minimum(0,h0)
h=h0+eta0
#h0[gnd>10]=gnd[gnd>10]
print('bat')
print(bat)
#plt.figure()
#plt.plot(-h0,'k')
#plt.plot(eta0,'y')
#plt.show()


x=1010
dx=10
dt=0.01
u=np.zeros((3,kcells+1))
eta=np.zeros((3,kcells))
eta[0]=eta0
g=9.81
ku=2*g*dt/dx
kh=2*dt/dx
eps=0.05
hmin=0.15
mu=dt*np.sqrt(g*h0)/dx
print('mu:', mu)
print('h0:', h0)

ini=0
x= kcells*dx
xi=np.arange(0,x,dx)
plt.figure()
#pasos de tiempo
for n in np.arange(85.5/dt):
    ni=int(n)%3
    #genera cond. ini.
    if n!=0:
        if n==1:
            u[1][1:-1]=u[0][1:-1]-ku*np.diff(eta[0])/2
            eta[1]=eta[0]-h*kh*np.diff(u[0])/2
        elif n>1:
            u[ni][1:-1]=u[ni-2][1:-1]-ku*np.diff(eta[ni-2])
            eta[ni]=eta[ni-2]-h*kh*np.diff(u[ni-2])
        u[ni][0]=0#u[ni][1]
        u[ni][-1]=0#u[ni][-2]
        idx_dry=np.asarray(h<hmin).nonzero()
        ucond=u[ni][1:]
        ucond[h<=hmin]=0
        #ucond[np.roll(h,-1)<=hmin]=0
        #ucond[np.roll(h,1)<=hmin]=0
        #print('u1:')
        #print(u[ni])
        u[ni][1:]=ucond
        #print('ucond:')
        #print(ucond)
        print('eta:',n*dt)
        print(eta[ni,0:50])
        print('ixd:',idx_dry)
        #print('where:')
        #print(np.where(ucond!=0))
        #print(np.where(h<=hmin))
            

        eta[ni][1:-1]=(1-eps)*eta[ni][1:-1]+0.5*eps*(np.roll(eta[ni],1)[1:-1]+np.roll(eta[ni],-1)[1:-1])
        eta[ni][idx_dry]=-h0[idx_dry]
        
        #eta[ni]-=np.minimum(0,h0)
        h=eta[ni]+h0
        #u[ni,0]=0
        #u[ni,-2]=0
        #u[ni,-1]=0

    if n%100==0:
        #wet=np.copy(h[ni])
        #wet[-bat>h]=np.NaN
        plt.plot(xi, -h0, 'b')
        plt.plot(xi, eta[ni],'.r')
        plt.axis([0, x, -1, 1.5])
        plt.title('t={}s, eps={}'.format(n*dt, eps))
        plt.plot(xi, u[ni][1:],'k')
        if n==0:
            plt.pause(1.5)
        else:
            plt.pause(0.3)
        plt.cla()
plt.show()

