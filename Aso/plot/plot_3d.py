#usage
#python plot_3d.py

from os import makedirs
from os.path import isdir,isfile

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

#figure
crater=[131.084935,32.884906,-1.506] #(lon,lat,dep)
lon_plt_min=130.95
lon_plt_max=131.15
lat_plt_min=32.8
lat_plt_max=33
lon_prfl=131.05
lat_prfl=32.9
cmap="rainbow_r"
cb_label="Vs (km/s)"
title="v01"

#velocity
vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/v01/vel_orig.npz")
idx=vel["idx"] #(lon,lat,dep)
idx_topo=vel["idx_topo"] #(lon,lat)
lon=vel["lon"]
lat=vel["lat"]
dep=vel["dep"]
res_inv=vel["res_inv"] #(lon,lat,dep)
num=vel["num"] #(lon,lat,dep)
vel=vel["vel"] #(lon,lat,dep)

idx_lon_prfl=round(res_inv[0]*(lon_prfl-lon[0]))
idx_lat_prfl=round(res_inv[1]*(lat_prfl-lat[0]))
vel[~idx]=np.nan
lat_lon=vel[np.arange(stop=num[0])[:,None],np.arange(stop=num[1]),idx_topo].T
lat_dep=vel[idx_lon_prfl]
dep_lon=vel[:,idx_lat_prfl].T
vmin=min(np.min(a=lat_lon),np.nanmin(a=lat_dep),np.nanmin(a=dep_lon))
vmax=max(np.max(a=lat_lon),np.nanmax(a=lat_dep),np.nanmax(a=dep_lon))
levels=np.linspace(start=vmin,stop=vmax,num=20)

#station
df=pd.read_csv(filepath_or_buffer="/home/tllc46/Aso/center,sep=" ",names=["stnm","stla","stlo","stel"])
df["stdp"]=-0.001*df["stel"]

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/figs/3d"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/v01.png"
if isfile(path=path_save):
    print("figure already exists")
    exit()

def init_topo():
    global lon_topo,topo_prfl_lon
    global lat_topo,topo_prfl_lat
    global topo

    topo=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/topo/Aso_topo_10m.npz")
    topo=topo["topo"] #(6000,6750)=(lat_topo,lon_topo) [m]
    lat_topo=np.arange(stop=6000)/9000+(32.5835-1/9000)
    lon_topo=np.arange(stop=6750)/9000+(130.7505-4/9000)

    idx_lon_min=int(np.floor(9000*(lon_plt_min-130.7505)+4))
    idx_lon_max=int(np.ceil(9000*(lon_plt_max-130.7505)+4))+1
    idx_lon_prfl=round(9000*(lon_prfl-130.7505)+4)
    lon_topo=lon_topo[idx_lon_min:idx_lon_max]

    idx_lat_min=int(np.floor(9000*(lat_plt_min-32.5835)+1))
    idx_lat_max=int(np.ceil(9000*(lat_plt_max-32.5835)+1))+1
    idx_lat_prfl=round(9000*(lat_prfl-32.5835)+1)
    lat_topo=lat_topo[idx_lat_min:idx_lat_max]

    topo_prfl_lon=-0.001*topo[idx_lat_prfl,idx_lon_min:idx_lon_max]
    topo_prfl_lat=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_prfl]
    topo=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_min:idx_lon_max]

def plot_3d():
    fig=plt.figure(figsize=(25.6,14.4))
    fig.set_layout_engine(layout="constrained")
    axes=fig.subplots(nrows=2,ncols=2,width_ratios=[2,1],height_ratios=[2,1])
    fig.delaxes(ax=axes[1,1])

    ax=axes[0,0]
    quadcontourset=ax.contourf(lon,lat,lat_lon,levels=levels,cmap=cmap)
    ax.contour(lon_topo,lat_topo,topo,colors="black",linestyles="solid")
    ax.axvline(x=lon_prfl,linestyle="dashed",color="black")
    ax.axhline(y=lat_prfl,linestyle="dashed",color="black")
    ax.scatter(x=df["stlo"],y=df["stla"],color="black",marker="^")
    ax.scatter(x=crater[0],y=crater[1],color="black",marker="*")

    ax.set_xlim(left=lon_plt_min,right=lon_plt_max)
    ax.set_ylim(bottom=lat_plt_min,top=lat_plt_max)
    ax.tick_params(labelbottom=False,labelleft=False)
    ax.grid(visible=True)

    ax=axes[0,1]
    ax.contourf(dep,lat,lat_dep,levels=levels,cmap=cmap)
    ax.axhline(y=lat_prfl,linestyle="dashed",color="black")
    ax.scatter(x=df["stdp"],y=df["stla"],color="black",marker="<")
    ax.scatter(x=crater[2],y=crater[1],color="black",marker="*")
    ax.plot(topo_prfl_lat,lat_topo,color="black")

    ax.set_xlabel(xlabel="depth (km)")
    ax.set_ylim(bottom=lat_plt_min,top=lat_plt_max)
    ax.set_ylabel(ylabel="latitude (°)")
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position(position="right")
    ax.grid(visible=True)

    ax=axes[1,0]
    ax.contourf(lon,dep,dep_lon,levels=levels,cmap=cmap)
    ax.axvline(x=lon_prfl,linestyle="dashed",color="black")
    ax.scatter(x=df["stlo"],y=df["stdp"],color="black",marker="^")
    ax.scatter(x=crater[0],y=crater[2],color="black",marker="*")
    ax.plot(lon_topo,topo_prfl_lon,color="black")

    ax.set_xlim(left=lon_plt_min,right=lon_plt_max)
    ax.set_xlabel(xlabel="longitude (°)")
    ax.set_ylabel(ylabel="depth (km)")
    ax.invert_yaxis()
    ax.grid(visible=True)

    ax=fig.add_axes(rect=(0.68,0.22,0.2,0.02))
    colorbar=fig.colorbar(mappable=quadcontourset,cax=ax,orientation="horizontal")
    colorbar.ax.xaxis.set_major_locator(locator=MaxNLocator(nbins=4))
    colorbar.set_label(label=cb_label)

    fig.suptitle(t=title)
    fig.savefig(fname=path_save)

def main():
    init_topo()
    plot_3d()

main()
