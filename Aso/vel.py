#usage
#v01_main.py → calc_tt.py → vel_reduce.py → calc_dtt.py

from os import makedirs
from os.path import isdir,isfile

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/v01"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/vel_orig.npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

topo=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/topo/Aso_topo_10m.npz")
topo=topo["topo"] #(6000,6750)=(lat_topo,lon_topo) [m]
lat_topo=np.arange(stop=6000)/9000+(32.5835-1/9000)
lon_topo=np.arange(stop=6750)/9000+(130.7505-4/9000)
interp_topo=RegularGridInterpolator(points=(lat_topo,lon_topo),values=topo)

lon_orig=np.array(object=[130.8,130.85,130.9,130.95,131,131.05,131.1,131.15,131.2,131.25,131.3,131.35,131.4])
lat_orig=np.array(object=[32.55,32.6,32.65,32.7,32.75,32.8,32.85,32.9,32.95,33,33.05,33.1,33.15])
ele_orig=np.array(object=[-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0])
vel_orig=np.empty(shape=(13,13,13)) #(lon_orig,lat_orig,ele_orig)

for i in range(13): #lon_orig
    for j in range(13): #lat_orig
        filepath="/home/tllc46/Documents/Aso/loc/velocity/"+f"{lat_orig[j]:.2f}-{lon_orig[i]:.2f}"
        df=pd.read_csv(filepath_or_buffer=filepath,sep="\\s+",names=["depth","velocity"])
        vel_orig[i,j]=np.flip(m=df["velocity"].to_numpy())
interp_vel=RegularGridInterpolator(points=(lon_orig,lat_orig,ele_orig),values=vel_orig)

lon_min=130.95
lon_max=131.15
lon_res_inv=2000
lon_res=1/lon_res_inv
lon_min_int=round(lon_min*lon_res_inv)
lon_max_int=round(lon_max*lon_res_inv)
lon=lon_res*np.arange(start=lon_min_int,stop=lon_max_int+1)
lon_num=lon_max_int+1-lon_min_int

lat_min=32.8
lat_max=33
lat_res_inv=2000
lat_res=1/lat_res_inv
lat_min_int=round(lat_min*lat_res_inv)
lat_max_int=round(lat_max*lat_res_inv)
lat=lat_res*np.arange(start=lat_min_int,stop=lat_max_int+1)
lat_num=lat_max_int+1-lat_min_int

ele_min=-2.6 #[km]
ele_max=1.6 #[km]
ele_res_inv=25
ele_res=1/ele_res_inv
ele_min_int=round(ele_min*ele_res_inv)
ele_max_int=round(ele_max*ele_res_inv)
ele=ele_res*np.arange(start=ele_min_int,stop=ele_max_int+1)
ele_num=ele_max_int+1-ele_min_int

lon_rdc_min=131.05
lon_rdc_max=131.1
idx_lon_min=round((lon_rdc_min-lon_min)*lon_res_inv)
idx_lon_max=round((lon_rdc_max-lon_min)*lon_res_inv)+1

lat_rdc_min=32.85
lat_rdc_max=32.925
idx_lat_min=round((lat_rdc_min-lat_min)*lat_res_inv)
idx_lat_max=round((lat_rdc_max-lat_min)*lat_res_inv)+1

vel=np.zeros(shape=(lon_num,lat_num,ele_num)) #(lon,lat,ele)
idx=np.zeros(shape=(lon_num,lat_num,ele_num),dtype=bool) #(lon,lat,ele)
idx_topo=np.empty(shape=(lon_num,lat_num),dtype=int) #(lon,lat)

for i in range(lon_num): #lon
    print(f"{i+1}/{lon_num}")
    for j in range(lat_num): #lat
        topo_ij=interp_topo(xi=(lat[j],lon[i]))
        topo_ij=ele_res*round(topo_ij*ele_res_inv/1000) #round to ele_res [km]
        for k in range(ele_num): #ele
            ele_k=ele[k] #[km]
            if topo_ij-6<=ele_k and ele_k<=topo_ij:
                vel[i,j,k]=interp_vel(xi=(lon[i],lat[j],ele_k-topo_ij))
                idx[i,j,k]=True
            elif topo_ij<ele_k:
                idx_topo[i,j]=ele_num-k
                break

dep=-np.flip(m=ele)
res_inv=np.array(object=[lon_res_inv,lat_res_inv,ele_res_inv])
res=np.array(object=[lon_res,lat_res,ele_res])
num=np.array(object=[lon_num,lat_num,ele_num])
idx_lon=np.array(object=[idx_lon_min,idx_lon_max])
idx_lat=np.array(object=[idx_lat_min,idx_lat_max])

vel[:,:,:]=np.flip(m=vel,axis=2) #(lon,lat,dep)
idx[:,:,:]=np.flip(m=idx,axis=2) #(lon,lat,dep)
np.savez(file=path_save,vel=vel,idx=idx,idx_topo=idx_topo,lon=lon,lat=lat,dep=dep,res=res,res_inv=res_inv,num=num,idx_lon=idx_lon,idx_lat=idx_lat)
