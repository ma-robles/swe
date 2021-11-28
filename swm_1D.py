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
h0=np.ones(bat.shape)*10
eta0=np.zeros(h0.shape)
#eta0[0:20]=1
eta0[40:61]=1
h=h0+eta0
#h0[gnd>10]=gnd[gnd>10]
print('bat')
print(bat)
#plt.figure()
#plt.plot(-h0,'k')
#plt.plot(eta0,'y')
#plt.show()

hmin=0.05

x=1010
dx=10
dt=0.01
u=np.zeros((3,kcells+1))
h=np.zeros((3,kcells))
eta=np.zeros((3,kcells))
eta[0]=eta0
h[0]=h0+eta0
g=9.81
ku=2*g*dt/dx
kh=2*dt/dx
eps=0.05
mu=dt*np.sqrt(g*h0)/dx
print('mu:', mu)
print('h0:', h0)

ini=0
x= kcells*dx
xi=np.arange(0,x,dx)
plt.figure()
f1= plt.subplot(2,1,1)
f2= plt.subplot(2,1,2)
#pasos de tiempo
for n in np.arange(24/dt):
    ni=int(n)%3
    #genera cond. ini.
    if n!=0:

        #u_anterior=u[ni-1]
        #u_siguiente=np.roll(u[ni-1],-1)
        #eta_anterior=np.roll(eta[ni-1],1)
        #eta_siguiente=eta[ni-1]
        #h_anterior=np.roll(h[ni-1],1)
        #h_siguiente=h[ni-1]
        if n==1:
            u[1][1:-1]=u[0][1:-1]-ku*np.diff(h[0])/2
            u[1][0]=u[1][1]
            u[1][-1]=u[1][-2]
            ucond=u[ni][1:]
            ucond[h[ni-1]<=hmin]=0
            print('ucond')
            print(ucond)
            #ucond[np.roll(h,-1)<=hmin]=0
            print(ucond)
            u[ni][1:]=ucond
            h[1]=h[0]-h[0]*kh*np.diff(u[0])/2
            print('u:', u[1])
            print('h:', h)
            print('eta:',eta[0],eta[1])
        elif n>1:
            u[ni][1:-1]=u[ni-2][1:-1]-ku*np.diff(h[ni-2])#(eta[ni-2])
            u[ni][0]=u[ni][1]
            u[ni][-1]=u[ni][-2]
            ucond=u[ni][1:]
            ucond[h[ni-1]<=hmin]=0
            #ucond[np.roll(h,-1)<=hmin]=0
            #ucond[np.roll(h,1)<=hmin]=0
            #print('u1:')
            #print(u[ni])
            u[ni][1:]=ucond
            #print('ucond:')
            #print(ucond)
            print('u:')
            print(u[ni,0:50])
            print('eta:',n*dt)
            print(h[ni,0:50])
            #print('where:')
            #print(np.where(ucond!=0))
            #print(np.where(h<=hmin))
            
            h[ni]=h[ni-2]-h[ni-2]*kh*np.diff(u[ni-2])

        h[ni][1:-1]=(1-eps)*h[ni][1:-1]+0.5*eps*(np.roll(h[ni],1)[1:-1]+np.roll(h[ni],-1)[1:-1])
        #eta[ni][1:-1]=(1-eps)*eta[ni][1:-1]+0.5*eps*(np.roll(eta[ni],1)[1:-1]+np.roll(eta[ni],-1)[1:-1])
        #h[ni]=(1-eps)*h[ni]+0.5*eps*(h_anterior+np.roll(h[ni-1],-1))
        #h=eta[ni]+h0
        #u[ni,0]=0
        #u[ni,-2]=0
        #u[ni,-1]=0

    if n%50==0:
        #wet=np.copy(h[ni])
        #wet[-bat>h]=np.NaN
        f1.plot(xi, -h0, 'b')
        f1.plot(xi, h[ni]-h0,'.r')
        f1.axis([0, x, -1, 1.5])
        f1.set_title('t={}s, eps={}'.format(n*dt, eps))
        f2.plot(xi, u[ni][1:])
        f2.axis([0, x, -1.0, 1.0])
        if n==0:
            plt.pause(1.5)
        else:
            plt.pause(0.3)
        f1.cla()
        f2.cla()
plt.show()

