from netCDF4 import Dataset
from matplotlib import pyplot as plt
import numpy as np

ifilename="isla_out.nc"

with Dataset(ifilename,"r") as ifile:
    h0=ifile["h0"][:]
    xi=ifile["x"][:]
    time=ifile["time"][:]
    plt.figure(figsize=(9,4))
    for i,t in enumerate(time):

        #crea η mas adecuada para graficación
        plt.plot(xi, -h0, 'b')
        plt.plot(xi, ifile["η"][i,:],'.r')
        plt.axis([0, xi[-1], -1, 1.5])
        plt.title('t={}s'.format(t))
        plt.plot(xi, ifile["u"][i,1:],'k')
        if i==0:
            plt.pause(1.5)
        else:
            plt.pause(0.1)
        plt.cla()
    plt.show()
