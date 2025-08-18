#usage
#v01_main.py → calc_tt.py → vel_reduce.py → calc_dtt.py
#python calc_tt.py v01 8

import sys
from os.path import isfile

import numpy as np
import pandas as pd

import pykonal

df=pd.read_csv(filepath_or_buffer="center",sep=" ",names=["stnm","stla","stlo","stel"])
df["stdp"]=-0.001*df["stel"] #station depth in [km]

idx_sta=int(sys.argv[2])
sta_coordinates=[df.loc[idx_sta,"stlo"],df.loc[idx_sta,"stla"],df.loc[idx_sta,"stdp"]]
stnm=df.loc[idx_sta,"stnm"]
path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/"+stnm+".npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/vel_orig.npz")
lon=vel["lon"]
lat=vel["lat"]
dep=vel["dep"]
res=vel["res"] #(lon,lat,dep)
num=vel["num"] #(lon,lat,dep)
idx_lon=vel["idx_lon"]
idx_lat=vel["idx_lat"]
vel=vel["vel"] #(lon,lat,dep)

travel_times=np.empty(shape=num) #(lon,lat,dep)

# Origin of the grid in spherical coordinates
origin_geo=[lat[-1],lon[0],dep[-1]]
origin_sph=pykonal.transformations.geo2sph(nodes=origin_geo) #radius (6371km earth radius),polar angle,azimuthal angle
r_min,theta_min,phi_min=origin_sph #r,θ,φ

# Extract resolutions
lon_res,lat_res,dep_res=res
lon_res=np.deg2rad(lon_res)
lat_res=np.deg2rad(lat_res)

# Initialize the Eikonal solver
solver=pykonal.solver.PointSourceSolver(coord_sys="spherical")

# Define the computational grid
solver.velocity.min_coords=origin_sph #r,θ,φ
solver.velocity.node_intervals=(dep_res,lat_res,lon_res) #dr,dθ,dφ
num_lon,num_lat,num_dep=num
solver.velocity.npts=(num_dep,num_lat,num_lon) #r,θ,φ
solver.velocity.values=np.flip(m=np.flip(m=np.swapaxes(a=vel,axis1=0,axis2=2),axis=1),axis=0) #(decreasing dep,decreasing lat,increasing lon)

# Convert the geographical coordinates of the receiver to spherical coordinates
sta_lon,sta_lat,sta_dep=sta_coordinates
sta_sph=pykonal.transformations.geo2sph(nodes=[sta_lat,sta_lon,sta_dep]) #r,θ,φ
r_sta,theta_sta,phi_sta=sta_sph #r,θ,φ
solver.src_loc=sta_sph #r,θ,φ

# Check if source is within grid:
if r_sta<r_min \
    or r_min+num_dep*dep_res<r_sta \
    or theta_sta<theta_min \
    or theta_min+num_lat*lat_res<theta_sta \
    or phi_sta<phi_min \
    or phi_min+num_lon*lon_res<phi_sta:
    print("Receiver is outside the grid! Pykonal will return Infs.")
    exit()

success=solver.solve()
if not success:
    print("solver execution went wrong")
    exit()

travel_times[:,:,:]=np.flip(m=np.flip(m=np.swapaxes(a=solver.traveltime.values,axis1=0,axis2=2),axis=2),axis=1) #(increasing lon,increasing lat,increasing dep)
travel_times=travel_times[idx_lon[0]:idx_lon[1],idx_lat[0]:idx_lat[1]]
travel_times[np.isinf(travel_times)]=0
print(stnm,np.sum(a=travel_times))
np.savez(file=path_save,travel_times=travel_times.flatten())
