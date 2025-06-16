#usage
#python tt.py 8

import sys

import numpy as np
import pandas as pd

import pykonal

df=pd.read_csv(filepath_or_buffer="center",sep=" ",names=["stnm","stla","stlo","stel"])
df["stdp"]=-0.001*df["stel"] #station depth in [km]

idx_sta=int(sys.argv[1])
sta_coordinates=[df.loc[idx_sta,"stlo"],df.loc[idx_sta,"stla"],df.loc[idx_sta,"stdp"]]

lat=0.002*np.arange(stop=101)+32.8
lon=0.002*np.arange(stop=101)+130.95
dep=0.1*np.arange(start=-16,stop=59) #[km]
vel=np.load(file="vel.npz")
vel=vel["vel"] #(lon,lat,dep)

travel_times=np.empty(shape=(101,101,75)) #(lon,lat,dep)

# Origin of the grid in spherical coordinates
origin_geo=[lat[-1],lon[0],dep[-1]]
origin_sph=pykonal.transformations.geo2sph(nodes=origin_geo) #radius (6371km earth radius),polar angle,azimuthal angle
r_min,theta_min,phi_min=origin_sph #r,θ,φ

# Extract resolutions
lon_res,lat_res,dep_res=(0.002,0.002,0.1)
lon_res=np.deg2rad(lon_res)
lat_res=np.deg2rad(lat_res)

# Initialize the Eikonal solver
solver=pykonal.solver.PointSourceSolver(coord_sys="spherical")

# Define the computational grid
solver.velocity.min_coords=origin_sph #r,θ,φ
solver.velocity.node_intervals=(dep_res,lat_res,lon_res) #dr,dθ,dφ
num_lon,num_lat,num_dep=(101,101,75)
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
travel_times[np.isinf(travel_times)]=0
np.savez(file=f"tt_{idx_sta:02}.npz",travel_times_flat=travel_times.flatten())
