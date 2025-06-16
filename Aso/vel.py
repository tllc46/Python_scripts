import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

topo=np.load(file="/home/tllc46/Downloads/topo/Aso_topo_10m.npz")
topo=topo["topo"] #(6000,6750)=(lat_topo,lon_topo) [m]
lat_topo=np.arange(stop=6000)/9000+(32.5835-1/9000)
lon_topo=np.arange(stop=6750)/9000+(130.7505-4/9000)
interp_topo=RegularGridInterpolator(points=(lat_topo,lon_topo),values=topo)

lat_vel=0.05*np.arange(stop=5)+32.8
lon_vel=0.05*np.arange(stop=5)+130.95
ele_vel=0.5*np.arange(start=-12,stop=1) #[km]
vel=np.empty(shape=(5,5,13)) #(lat_vel,lon_vel,ele_vel)

for i in range(5): #lat_vel
    for j in range(5): #lon_vel
        filepath="/home/tllc46/Documents/Aso/loc/velocity/"+f"{lat_vel[i]:.2f}-{lon_vel[j]:.2f}"
        df=pd.read_csv(filepath_or_buffer=filepath,sep="\\s+",names=["depth","velocity"])
        vel[i,j]=np.flip(m=df["velocity"].to_numpy())
interp_vel=RegularGridInterpolator(points=(lat_vel,lon_vel,ele_vel),values=vel)

lat_grid=0.002*np.arange(stop=101)+32.8
lon_grid=0.002*np.arange(stop=101)+130.95
ele_grid=0.1*np.arange(start=-58,stop=17) #[km]
vel_grid=np.zeros(shape=(101,101,75)) #(lat_grid,lon_grid,ele_grid)
idx_grid=np.zeros(shape=(101,101,75),dtype=bool) #(lat_grid,lon_grid,ele_grid)

for i in range(101): #lat_grid
    for j in range(101): #lon_grid
        topo_grid=interp_topo(xi=(lat_grid[i],lon_grid[j]))
        topo_grid=0.1*round(topo_grid/100) #round to ele_grid resolution [km]
        for k in range(75): #ele_grid
            ele_k=ele_grid[k] #[km]
            if topo_grid-6<=ele_k and ele_k<=topo_grid:
                vel_grid[i,j,k]=interp_vel(xi=(lat_grid[i],lon_grid[j],ele_k-topo_grid))
                idx_grid[i,j,k]=True

vel_grid=np.swapaxes(a=vel_grid,axis1=0,axis2=1) #(lon_grid,lat_grid,ele_grid)
vel_grid=np.flip(m=vel_grid,axis=2) #(lon_grid,lat_grid,dep_grid)
idx_grid=np.swapaxes(a=idx_grid,axis1=0,axis2=1) #(lon_grid,lat_grid,ele_grid)
idx_grid=np.flip(m=idx_grid,axis=2) #(lon_grid,lat_grid,dep_grid)
np.savez(file="vel.npz",vel=vel_grid,idx_flat=idx_grid.flatten())
