#usage
#v06_main.py → set_grid.py → calc_tt.py → calc_dtt.py

from os import makedirs
from os.path import isdir,isfile

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

from aso98_extend import lon_orig,lat_orig,ele_orig,vel_orig

vel_orig/=1.704

path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/v06"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/vel.npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

topo=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/topo/Aso_topo_10m.npz")
topo=topo["topo"] #(6000,6750)=(lat_topo,lon_topo) [m]
lat_topo=np.arange(stop=6000)/9000+(32.5835-1/9000)
lon_topo=np.arange(stop=6750)/9000+(130.7505-4/9000)
interp_topo=RegularGridInterpolator(points=(lat_topo,lon_topo),values=topo)

interp_vel=RegularGridInterpolator(points=(lon_orig,lat_orig,ele_orig),values=vel_orig,bounds_error=False,fill_value=None)

lon_min=131.05
lon_max=131.15
lon_res_inv=10000
lon_res=1/lon_res_inv
lon_min_int=round(lon_min*lon_res_inv)
lon_max_int=round(lon_max*lon_res_inv)+1
lon=lon_res*np.arange(start=lon_min_int,stop=lon_max_int)
lon_num=lon_max_int-lon_min_int

lat_min=32.84
lat_max=32.93
lat_res_inv=10000
lat_res=1/lat_res_inv
lat_min_int=round(lat_min*lat_res_inv)
lat_max_int=round(lat_max*lat_res_inv)+1
lat=lat_res*np.arange(start=lat_min_int,stop=lat_max_int)
lat_num=lat_max_int-lat_min_int

ele_min=0 #[km]
ele_max=1.6 #[km]
ele_res_inv=100
ele_res=1/ele_res_inv
ele_min_int=round(ele_min*ele_res_inv)
ele_max_int=round(ele_max*ele_res_inv)+1
ele=ele_res*np.arange(start=ele_min_int,stop=ele_max_int)
ele_num=ele_max_int-ele_min_int

vel=np.zeros(shape=(lon_num,lat_num,ele_num)) #(lon,lat,ele)
idx=np.zeros(shape=(lon_num,lat_num,ele_num),dtype=bool) #(lon,lat,ele)

for i in range(lon_num): #lon
    print(f"{i+1}/{lon_num}")
    for j in range(lat_num): #lat
        topo_ij=interp_topo(xi=(lat[j],lon[i]))
        topo_ij=ele_res*round(topo_ij*ele_res_inv/1000) #round to ele_res [km]
        for k in range(ele_num): #ele
            vel[i,j,k]=interp_vel(xi=(lon[i],lat[j],ele[k]))
            if ele[k]<=topo_ij:
                idx[i,j,k]=True

dep=-np.flip(m=ele)
res_inv=np.array(object=[lon_res_inv,lat_res_inv,ele_res_inv])
res=np.array(object=[lon_res,lat_res,ele_res])
num=np.array(object=[lon_num,lat_num,ele_num])

vel[:,:,:]=np.flip(m=vel,axis=2) #(lon,lat,dep)
idx[:,:,:]=np.flip(m=idx,axis=2) #(lon,lat,dep)
np.savez(file=path_save,vel=vel,idx=idx,lon=lon,lat=lat,dep=dep,res=res,res_inv=res_inv,num=num)
