#usage
#v01_main.py → set_grid.py → calc_tt.py → calc_dtt.py
#python set_grid.py v06

import sys
from os import makedirs
from os.path import isdir,isfile

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/grid.npz"
#if isfile(path=path_save):
#    print("data already exists")
#    exit()

topo=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/topo/Aso_topo_10m.npz")
topo=topo["topo"] #(6000,6750)=(lat_topo,lon_topo) [m]
lat_topo=np.arange(stop=6000)/9000+(32.5835-1/9000)
lon_topo=np.arange(stop=6750)/9000+(130.7505-4/9000)
interp_topo=RegularGridInterpolator(points=(lat_topo,lon_topo),values=topo)

vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/vel.npz")
lon=vel["lon"]
lat=vel["lat"]
dep=vel["dep"]
res_inv=vel["res_inv"]

lon_min=131.06
lon_max=131.15
lon_res_inv=2000
lon_res=1/lon_res_inv
idx_lon_min=round((lon_min-lon[0])*res_inv[0])
idx_lon_max=round((lon_max-lon[0])*res_inv[0])+1
lon_step=res_inv[0]//lon_res_inv

lon_min_int=round(lon_min*lon_res_inv)
lon_max_int=round(lon_max*lon_res_inv)+1
lon=lon_res*np.arange(start=lon_min_int,stop=lon_max_int)
lon_num=lon_max_int-lon_min_int

lat_min=32.85
lat_max=32.91
lat_res_inv=2000
lat_res=1/lat_res_inv
idx_lat_min=round((lat_min-lat[0])*res_inv[1])
idx_lat_max=round((lat_max-lat[0])*res_inv[1])+1
lat_step=res_inv[1]//lat_res_inv

lat_min_int=round(lat_min*lat_res_inv)
lat_max_int=round(lat_max*lat_res_inv)+1
lat=lat_res*np.arange(start=lat_min_int,stop=lat_max_int)
lat_num=lat_max_int-lat_min_int

dep_min=-1.6
dep_max=0
dep_res_inv=25
dep_res=1/dep_res_inv
idx_dep_min=round((dep_min-dep[0])*res_inv[2])
idx_dep_max=round((dep_max-dep[0])*res_inv[2])+1
dep_step=res_inv[2]//dep_res_inv

dep_min_int=round(dep_min*dep_res_inv)
dep_max_int=round(dep_max*dep_res_inv)+1
dep=dep_res*np.arange(start=dep_min_int,stop=dep_max_int)
dep_num=dep_max_int-dep_min_int

ele_res_inv=dep_res_inv
ele_res=1/ele_res_inv
ele=-np.flip(m=dep)
ele_num=dep_num

idx=np.zeros(shape=(lon_num,lat_num,ele_num),dtype=bool) #(lon,lat,ele)
idx_topo=np.zeros(shape=(lon_num,lat_num),dtype=int) #(lon,lat)

for i in range(lon_num): #lon
    for j in range(lat_num): #lat
        topo_ij=interp_topo(xi=(lat[j],lon[i]))
        topo_ij=ele_res*round(topo_ij*ele_res_inv/1000) #round to ele_res [km]
        for k in range(ele_num): #ele
            if ele[k]<=topo_ij:
                idx[i,j,k]=True
            else:
                idx_topo[i,j]=ele_num-k
                break

res_inv=np.array(object=[lon_res_inv,lat_res_inv,ele_res_inv])
res=np.array(object=[lon_res,lat_res,ele_res])
num=np.array(object=[lon_num,lat_num,ele_num])
idx_lon=np.array(object=[idx_lon_min,idx_lon_max,lon_step])
idx_lat=np.array(object=[idx_lat_min,idx_lat_max,lat_step])
idx_dep=np.array(object=[idx_dep_min,idx_dep_max,dep_step])

idx[:,:,:]=np.flip(m=idx,axis=2) #(lon,lat,dep)
np.savez(file=path_save,idx=idx,idx_topo=idx_topo,lon=lon,lat=lat,dep=dep,res=res,res_inv=res_inv,num=num,idx_lon=idx_lon,idx_lat=idx_lat,idx_dep=idx_dep)
