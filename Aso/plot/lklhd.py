#usage
#python lklhd.py t05 x195 v06 b03 2015-01-15T03:00:00

import sys
from os import makedirs
from os.path import isdir,isfile
from importlib import import_module
from datetime import datetime

import numpy as np
from scipy.signal import butter,sosfilt,hilbert
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

import fncts_topo_10 as ft

bool_fltr=False

sys.path.append("/home/tllc46/Aso/tstep")
sys.path.append("/home/tllc46/Aso/xcorr/params")
sys.path.append("/home/tllc46/Aso/bp/params")

mdl_t=import_module(name=sys.argv[1])
mdl_x=import_module(name=sys.argv[2])
mdl_b=import_module(name=sys.argv[4])

#sampling rate
sampling_rate=100 #[Hz]

#average window
dt=datetime.strptime(sys.argv[5],"%Y-%m-%dT%H:%M:%S")
mdl_t.parse_dt(dt=dt)

#velocity grid
grid=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/grid.npz")
idx=grid["idx"] #(lon,lat,dep)
idx_topo=grid["idx_topo"] #(lon,lat)
lon=grid["lon"]
lat=grid["lat"]
dep=grid["dep"]
num=grid["num"] #(lon,lat,dep)
idx_flat=idx.flatten() #(nnode,)
nnode=np.prod(a=num)

#cross correlation
xcorr=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/xcorr/"+sys.argv[2]+"."+mdl_t.name+".npz")
xcorr=xcorr["xcorr"][mdl_t.idx_avg] #(ntriu,npts_sub)
npts_sub=xcorr.shape[1]

#travel time difference
idx_dtt=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/idx_dtt.npz")
idx_dtt=idx_dtt["idx_dtt"] #(ntriu,nnode)
idx_dtt+=npts_sub//2

#station
df=pd.read_csv(filepath_or_buffer=mdl_x.info_sta,sep=" ",names=["stnm","stla","stlo","stel"])
df["stdp"]=-0.001*df["stel"]
nsta=len(df)
idx_triu=np.triu_indices(n=nsta,k=1)
stnm_pairs=[]
for i in range(nsta):
    for j in range(i+1,nsta):
        stnm_pairs.append(df.loc[i,"stnm"]+"-"+df.loc[j,"stnm"])
ntriu=nsta*(nsta-1)//2

#station pair choice
choice_single=np.ones(shape=nsta,dtype=bool)
for i in mdl_b.exclude_single:
    idx_stnm=df[df["stnm"]==i].index[0]
    choice_single[idx_stnm]=False
choice_pair=(choice_single[:,None] & choice_single)[idx_triu] #(ntriu,)
for i in mdl_b.exclude_pair:
    idx_pair=stnm_pairs.index(i)
    choice_pair[idx_pair]=False

#bandpass filter
if bool_fltr:
    sos_butter=butter(N=4,Wn=[3,6],btype="band",fs=sampling_rate,output="sos")

#array
xcorr_pair=np.empty(shape=npts_sub)
beam=np.zeros(shape=nnode)
idx_loc=np.empty(shape=3,dtype=int)

#figure
crater=[131.085,32.885] #(lon,lat)
coord_min=[131.05,32.84]
coord_max=[131.15,32.93]
idx_min=[123,120,23] #(lon,lat,dep)

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/figs/beam"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+dt.strftime(format="%Y-%m-%dT%H%M%S")+".png"

for i in range(ntriu):
    xcorr_pair[:]=xcorr[i]
    if np.isnan(xcorr_pair[0]) or not choice_pair[i]:
        continue

    if bool_fltr:
        xcorr_pair[:]=sosfilt(sos=sos_butter,x=xcorr_pair)
        xcorr_pair[:]=sosfilt(sos=sos_butter,x=xcorr_pair[::-1])[::-1]
    xcorr_pair[:]=abs(hilbert(x=xcorr_pair))
    xcorr_pair[:]=gaussian_filter1d(input=xcorr_pair,sigma=mdl_b.sigma)
    if mdl_b.normalize:
        xcorr_pair/=max(xcorr_pair)

    beam+=xcorr_pair[idx_dtt[i]]

beam_full=np.copy(a=beam)
beam[~idx_flat]=0
idx_max=np.argmax(a=beam)
beam_full/=beam_full[idx_max]
beam_full=np.reshape(beam_full,shape=num)
idx_loc[:]=np.unravel_index(indices=idx_max,shape=num)

lat_lon=beam_full[np.arange(stop=num[0])[:,None],np.arange(stop=num[1]),idx_topo].T
lat_dep=np.copy(a=beam_full[idx_loc[0]])
dep_lon=np.copy(a=beam_full[:,idx_loc[1]])
lat_dep[~idx[idx_loc[0]]]=np.nan
dep_lon[~idx[:,idx_loc[1]]]=np.nan
vmin=min(np.min(a=lat_lon),np.nanmin(a=lat_dep),np.nanmin(a=dep_lon))
lat_dep=beam_full[idx_loc[0]]
dep_lon=beam_full[:,idx_loc[1]].T

ft.init_lklhd(coord_min=coord_min,coord_max=coord_max,idx_loc=idx_loc[:2],grid=grid)

levels=np.linspace(start=vmin,stop=1,num=20)

#fig=plt.figure(figsize=(8.27,4.65),layout="constrained")
fig=plt.figure(figsize=(25.6,14.4),layout="constrained")
axes=fig.subplots(nrows=2,ncols=2,sharex="col",sharey="row",width_ratios=[2,1],height_ratios=[2,1])
fig.delaxes(ax=axes[1,1])

ax=axes[0,0]
quadcontourset=ax.contourf(lon,lat,lat_lon,levels=levels,cmap="rainbow")
ax.contour(ft.lon,ft.lat,ft.topo,colors="black",linestyles="solid")
ax.axvline(x=lon[idx_loc[0]],linestyle="dashed",color="black")
ax.axhline(y=lat[idx_loc[1]],linestyle="dashed",color="black")
ax.scatter(x=df["stlo"],y=df["stla"],color="black",marker="^")
ax.scatter(x=crater[0],y=crater[1],color="black",marker="*")
ax.scatter(x=lon[idx_loc[0]],y=lat[idx_loc[1]],color="None",edgecolor="black")
ax.scatter(x=lon[idx_min[0]],y=lat[idx_min[1]],color="None",marker="s",edgecolor="black")

ax.set_xlim(left=coord_min[0],right=coord_max[0])
ax.set_ylim(bottom=coord_min[1],top=coord_max[1])
ax.tick_params(labelleft=False)
ax.grid(visible=True)

ax=axes[0,1]
ax.contourf(dep,lat,lat_dep,levels=levels,cmap="rainbow")
ax.fill_betweenx(y=ft.lat,x1=dep[0],x2=ft.topo_prfl_lat,color="white")
ax.axhline(y=lat[idx_loc[1]],linestyle="dashed",color="black")
ax.scatter(x=df["stdp"],y=df["stla"],color="black",marker="<")
ax.plot(ft.topo_prfl_lat,ft.lat,color="black")
ax.scatter(x=dep[idx_loc[2]],y=lat[idx_loc[1]],color="None",edgecolor="black")
ax.scatter(x=dep[idx_min[2]],y=lat[idx_min[1]],color="None",marker="s",edgecolor="black")

ax.set_xlim(left=dep[0],right=dep[-1])
ax.set_xlabel(xlabel="depth (km)")
ax.set_ylabel(ylabel="latitude (°)")
ax.yaxis.tick_right()
ax.yaxis.set_label_position(position="right")
ax.tick_params(labelbottom=True,labelright=True)
ax.grid(visible=True)

ax=axes[1,0]
ax.contourf(lon,dep,dep_lon,levels=levels,cmap="rainbow")
ax.fill_between(x=ft.lon,y1=dep[0],y2=ft.topo_prfl_lon,color="white")
ax.axvline(x=lon[idx_loc[0]],linestyle="dashed",color="black")
ax.scatter(x=df["stlo"],y=df["stdp"],color="black",marker="^")
ax.plot(ft.lon,ft.topo_prfl_lon,color="black")
ax.scatter(x=lon[idx_loc[0]],y=dep[idx_loc[2]],color="None",edgecolor="black")
ax.scatter(x=lon[idx_min[0]],y=dep[idx_min[2]],color="None",marker="s",edgecolor="black")

ax.set_xlabel(xlabel="longitude (°)")
ax.set_ylim(bottom=dep[0],top=dep[-1])
ax.set_ylabel(ylabel="depth (km)")
ax.invert_yaxis()
ax.grid(visible=True)

ax=fig.add_axes(rect=(0.68,0.22,0.2,0.02))
colorbar=fig.colorbar(mappable=quadcontourset,cax=ax,orientation="horizontal")
colorbar.ax.xaxis.set_major_locator(locator=MaxNLocator(nbins=4))
colorbar.set_label(label="beam")

fig.savefig(fname=path_save)
