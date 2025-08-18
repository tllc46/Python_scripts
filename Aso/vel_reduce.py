#usage
#v01_main.py → calc_tt.py → vel_reduce.py → calc_dtt.py
#python vel_reduce.py v01

import sys
from os.path import isfile

import numpy as np

path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/vel.npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/vel_orig.npz")
idx=vel["idx"]
idx_topo=vel["idx_topo"]
lon=vel["lon"]
lat=vel["lat"]
dep=vel["dep"]
res=vel["res"]
res_inv=vel["res_inv"]
num=vel["num"]
idx_lon=vel["idx_lon"]
idx_lat=vel["idx_lat"]
vel=vel["vel"]

idx=idx[idx_lon[0]:idx_lon[1],idx_lat[0]:idx_lat[1]]
idx_topo=idx_topo[idx_lon[0]:idx_lon[1],idx_lat[0]:idx_lat[1]]
lon=lon[idx_lon[0]:idx_lon[1]]
lat=lat[idx_lat[0]:idx_lat[1]]
vel=vel[idx_lon[0]:idx_lon[1],idx_lat[0]:idx_lat[1]]
num=np.array(object=[idx_lon[1]-idx_lon[0],idx_lat[1]-idx_lat[0],num[2]])

np.savez(file=path_save,vel=vel,idx=idx,idx_topo=idx_topo,lon=lon,lat=lat,dep=dep,res=res,res_inv=res_inv,num=num)
