import numpy as np
from netCDF4 import Dataset

ifilename="nc/ladera.nc"
ofilename="ladera_out.nc"

with Dataset(ifilename,"r") as ifile:
    η0=ifile["η0"][:]
    h0=ifile["h0"][:]
    xi=ifile["x"][:]
    yi=ifile["y"][:]

#tipo de frontera
#frontera para isla
boundary_open= False
#constante para filtro
ϵ=0.03
#umbral para celdas secas
#hmin para isla
hmin=0.15
dx=5
dy=5
dt=0.01
g=9.81
ku=2*g*dt/dx
kv=2*g*dt/dy
khx=2*dt/dx
khy=2*dt/dy
xcells= h0.shape[0]
ycells= h0.shape[1]
u=np.zeros((3,xcells+1))
v=np.zeros((3,ycells+1))
η=np.zeros((3, xcells, ycells))
η[0]=η0
h=h0+η0
#μ=dt*np.sqrt(g*np.max(h0))/dx
#print('μ=', μ)
#tiempo de simulación[s]
time_tot=400.5
#divisor de tiempo para almacenar
time_div=100

with Dataset(ofilename, 'w') as ofile:
    ofile.description="salida de script swe 2D"
    dim_t= ofile.createDimension("time", None)
    dim_x= ofile.createDimension("x", xcells)
    dim_y= ofile.createDimension("y", ycells)
    dim_xstag= ofile.createDimension("x_stag", xcells+1)
    dim_ystag= ofile.createDimension("y_stag", ycells+1)

    var_t= ofile.createVariable('time', 'f8', ("time",))
    var_x= ofile.createVariable('x', 'f8', ("x",))
    var_y= ofile.createVariable('y', 'f8', ("y",))
    var_h0= ofile.createVariable('h0', 'f8', ("x","y"))
    var_η= ofile.createVariable('η', 'f8', ("time","x", "y"))
    var_u= ofile.createVariable('u', 'f8', 
            ("time", "x_stag", "y"))
    var_v= ofile.createVariable('v', 'f8', 
            ("time", "x", "y_stag"))
    var_x[:]=xi
    var_y[:]=yi
    var_h0[:]=h0
    print(var_t.shape, var_η.shape, var_u.shape)
    for n in np.arange(time_tot/dt):
        ni=int(n)%3
        if n!=0:
            if n==1:
                u[1][1:-1]=u[0][1:-1]-ku*np.diff(η[0])/2
                hu=h*u[0][1:]
                hu=np.insert(hu,0, 0)
                η[1]=η[0]-kh*np.diff(hu)/2
            elif n>1:
                u[ni][1:-1]=u[ni-2][1:-1,:]-ku*np.diff(η[ni-2])
                hu=h*u[ni-2][1:,:]
                hu=np.insert(hu,0, np.zeros(ycells), axis=1)
                η[ni]=η[ni-2]-kh*np.diff(hu)
            if boundary_open==False:
                u[ni][0]=0
                u[ni][-1]=0
            else:
                u[ni][0]=u[ni][1]
                u[ni][-1]=u[ni][-2]
            h=η[ni]+h0
            #elimina cambios en u para celdas no húmedas
            ucond=u[ni][1:]
            ucond[h<=hmin]=0
            u[ni][1:]=ucond
                
            #filtra para suavizar η
            #vecinos
            η_v1=np.roll(η[ni],1)[1:-1]
            η_v2=np.roll(η[ni],-1)[1:-1]
            η[ni][1:-1]=(1-ϵ)*η[ni][1:-1]+ 0.5*ϵ*(η_v1+η_v2)
        #cálculo de celdas secas
        idx_dry=np.asarray(h<hmin).nonzero()[0]
        #celdas límites
        lim=(i for i in idx_dry if (h[i-1]>=hmin or h[i+1]>=hmin) )
        try:
            lim = list(lim)
        except IndexError:
            pass
        #celdas secas
        dry = ( i for i in idx_dry if i not in lim)
        dry = list(dry)
        #elimina cambios en η para celdas secas
        η[ni][dry]=-h0[dry]
        h=η[ni]+h0
        if int(n)%time_div==0:
            η_plot=np.copy(η[ni])
            η_plot[dry]=np.NaN
            n2=n//time_div
            var_t[n2]=n*dt
            var_η[n2,:]=η_plot
            #var_η[n2,:]=η[ni]
            var_u[n2,:]=u[ni]
            #print(n, var_t.shape, var_η.shape, var_u.shape)
