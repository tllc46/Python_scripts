#usage
#python vel_3d.py v06

import sys

from os import makedirs
from os.path import isdir,isfile

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import fncts_topo_10 as ft

bool_ray=False

#figure
crater=[131.085,32.885] #(lon,lat)
coord_min=[131.05,32.84]
coord_max=[131.15,32.93]
coord_prfl=crater

#velocity
vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/vel.npz")
idx=vel["idx"] #(lon,lat,dep)
lon=vel["lon"]
lat=vel["lat"]
dep=vel["dep"]
res_inv=vel["res_inv"] #(lon,lat,dep)
vel=vel["vel"]

idx_lon_prfl=round(res_inv[0]*(coord_prfl[0]-lon[0]))
idx_lat_prfl=round(res_inv[1]*(coord_prfl[1]-lat[0]))

lat_lon=vel[:,:,100].T
lat_dep=np.copy(a=vel[idx_lon_prfl])
dep_lon=np.copy(a=vel[:,idx_lat_prfl])
lat_dep[~idx[idx_lon_prfl]]=np.nan
dep_lon[~idx[:,idx_lat_prfl]]=np.nan
vmin=min(np.min(a=lat_lon),np.nanmin(a=lat_dep),np.nanmin(a=dep_lon))
vmax=max(np.max(a=lat_lon),np.nanmax(a=lat_dep),np.nanmax(a=dep_lon))
lat_dep=vel[idx_lon_prfl]
dep_lon=vel[:,idx_lat_prfl].T

#station
df=pd.read_csv(filepath_or_buffer="/home/tllc46/Aso/center",sep=" ",names=["stnm","stla","stlo","stel"])
df["stdp"]=-0.001*df["stel"]
nsta=len(df)

#ray
if bool_ray:
    rays=[]
    for i in range(nsta):
        ray=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/rays/"+df.loc[i,"stnm"]+".npz")
        rays.append(ray["ray"])

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/figs/3d"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[1]+".png"

ft.init_prfl(coord_min=coord_min,coord_max=coord_max,coord_prfl=coord_prfl)

#fig=plt.figure(figsize=(8.27,4.65),layout="constrained")
fig=plt.figure(figsize=(25.6,14.4),layout="constrained")
axes=fig.subplots(nrows=2,ncols=2,sharex="col",sharey="row",width_ratios=[2,1],height_ratios=[2,1])
fig.delaxes(ax=axes[1,1])

ax=axes[0,0]
quadmesh=ax.pcolormesh(lon,lat,lat_lon,cmap="rainbow_r",vmin=vmin,vmax=vmax)
ax.contour(ft.lon,ft.lat,ft.topo,colors="black",linestyles="solid")
if bool_ray:
    for ray in rays:
        ax.plot(ray[:,1],ray[:,0],color="black")
ax.axvline(x=coord_prfl[0],linestyle="dashed",color="black")
ax.axhline(y=coord_prfl[1],linestyle="dashed",color="black")
ax.scatter(x=df["stlo"],y=df["stla"],color="black",marker="^")
ax.scatter(x=crater[0],y=crater[1],color="black",marker="*")

ax.set_xlim(left=coord_min[0],right=coord_max[0])
ax.set_ylim(bottom=coord_min[1],top=coord_max[1])
ax.tick_params(labelleft=False)
ax.grid(visible=True)

ax=axes[0,1]
ax.pcolormesh(dep,lat,lat_dep,cmap="rainbow_r",vmin=vmin,vmax=vmax)
ax.fill_betweenx(y=ft.lat,x1=dep[0],x2=ft.topo_prfl_lat,color="white")
if bool_ray:
    for ray in rays:
        ax.plot(ray[:,2],ray[:,0],color="black")
ax.axhline(y=coord_prfl[1],linestyle="dashed",color="black")
ax.axvline(x=-0.6,linestyle="dashed",color="black")
ax.scatter(x=df["stdp"],y=df["stla"],color="black",marker="<")
ax.plot(ft.topo_prfl_lat,ft.lat,color="black")

ax.set_xlim(left=dep[0],right=dep[-1])
ax.set_xlabel(xlabel="depth (km)")
ax.set_ylabel(ylabel="latitude (°)")
ax.yaxis.tick_right()
ax.yaxis.set_label_position(position="right")
ax.tick_params(labelbottom=True,labelright=True)
ax.grid(visible=True)

ax=axes[1,0]
ax.pcolormesh(lon,dep,dep_lon,cmap="rainbow_r",vmin=vmin,vmax=vmax)
ax.fill_between(x=ft.lon,y1=dep[0],y2=ft.topo_prfl_lon,color="white")
if bool_ray:
    for ray in rays:
        ax.plot(ray[:,1],ray[:,2],color="black")
ax.axvline(x=coord_prfl[0],linestyle="dashed",color="black")
ax.axhline(y=-0.6,linestyle="dashed",color="black")
ax.scatter(x=df["stlo"],y=df["stdp"],color="black",marker="^")
ax.plot(ft.lon,ft.topo_prfl_lon,color="black")

ax.set_xlabel(xlabel="longitude (°)")
ax.set_ylim(bottom=dep[0],top=dep[-1])
ax.set_ylabel(ylabel="depth (km)")
ax.invert_yaxis()
ax.grid(visible=True)

ax=fig.add_axes(rect=(0.68,0.22,0.2,0.02))
colorbar=fig.colorbar(mappable=quadmesh,cax=ax,orientation="horizontal")
colorbar.set_label(label="Vs (km/s)")

fig.savefig(fname=path_save)
