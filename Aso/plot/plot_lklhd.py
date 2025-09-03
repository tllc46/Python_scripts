#usage
#python plot_lklhd.py t01 x02 v01 b04 2016-10-19T03:00:00

import sys
from os import makedirs
from os.path import isdir,isfile
from importlib import import_module
from datetime import datetime,timedelta

import numpy as np
from scipy.signal import butter,sosfilt,hilbert
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

bool_fltr=False

sys.path.append("/home/tllc46/Aso/tstep")
sys.path.append("/home/tllc46/Aso/xcorr")
sys.path.append("/home/tllc46/Aso/bp")

mdl_t=import_module(name=sys.argv[1])
mdl_x=import_module(name=sys.argv[2])
mdl_b=import_module(name=sys.argv[4])

#constant
sec_day=86400 #[s]

#sampling rate
sampling_rate=100

#sub window
len_sub=48
npts_sub=len_sub*sampling_rate

#average window
shift_avg=mdl_t.shift_avg
offset_avg=mdl_t.offset_avg
dt=datetime.strptime(sys.argv[5],"%Y-%m-%dT%H:%M:%S")
dt_start_avg=datetime(year=dt.year,month=dt.month,day=dt.day)+timedelta(seconds=offset_avg)
if dt<dt_start_avg:
    print("time is before the average window offset")
    exit()
start_avg=(dt-dt_start_avg).seconds
if start_avg%shift_avg:
    print("time matching average window doesn't exist")
    exit()
navg_day=sec_day//shift_avg
mdl_t.parse_dt(dt=dt)
idx_avg=(dt-mdl_t.dt_name).days*navg_day+start_avg//shift_avg

#velocity
vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/vel.npz")
idx_flat=vel["idx"].flatten() #(lon,lat,dep) → (nnode,)
idx_topo=vel["idx_topo"] #(lon,lat)
lon=vel["lon"]
lat=vel["lat"]
dep=vel["dep"]
res_inv=vel["res_inv"] #(lon,lat,dep)
num=vel["num"] #(lon,lat,dep)
nnode=np.prod(a=num)

#travel time difference
idx_dtt=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/idx_dtt.npz")
idx_dtt=idx_dtt["idx_dtt"] #(ntriu,nnode)

#cross correlation
xcorr=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/xcorr/"+sys.argv[2]+"."+mdl_t.name+".npz")
xcorr=xcorr["xcorr"][idx_avg] #(ntriu,npts_sub)

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
choice_pair=(choice_single[:,None] & choice_single)[idx_triu]
for i in mdl_b.exclude_pair:
    idx_pair=stnm_pairs.index(i)
    choice_pair[idx_pair]=False

#bandpass filter
if bool_fltr:
    sos_butter=butter(N=4,Wn=[3,6],btype="band",fs=sampling_rate,output="sos")

#array
xcorr_pair=np.empty(shape=npts_sub)
beam=np.empty(shape=nnode)
idx_loc=np.empty(shape=3,dtype=int)

#figure
crater=[131.084935,32.884906,-1.506] #(lon,lat,dep)
lon_plt_min=130.95
lon_plt_max=131.15
lat_plt_min=32.8
lat_plt_max=33

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/figs/beam"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+dt.strftime(format="%Y-%m-%dT%H%M%S")+".png"
if isfile(path=path_save):
    print("figure already exists")
    exit()

def beamform():
    global xcorr_pair,beam
    nstack=0
    beam[:]=0

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
        nstack+=1

    beam[~idx_flat]=np.nan
    beam/=nstack
    idx_max=np.nanargmax(a=beam)
    beam=np.reshape(beam,shape=num)
    idx_loc[:]=np.unravel_index(indices=idx_max,shape=num)

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
    lon_inc=9000/res_inv[0]
    idx_lon_prfl=round(9000*(lon[0]-130.7505)+4+lon_inc*idx_loc[0])
    lon_topo=lon_topo[idx_lon_min:idx_lon_max]

    idx_lat_min=int(np.floor(9000*(lat_plt_min-32.5835)+1))
    idx_lat_max=int(np.ceil(9000*(lat_plt_max-32.5835)+1))+1
    lat_inc=9000/res_inv[1]
    idx_lat_prfl=round(9000*(lat[0]-32.5835)+1+lat_inc*idx_loc[1])
    lat_topo=lat_topo[idx_lat_min:idx_lat_max]

    topo_prfl_lon=-0.001*topo[idx_lat_prfl,idx_lon_min:idx_lon_max]
    topo_prfl_lat=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_prfl]
    topo=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_min:idx_lon_max]

def plot_3d():
    vmin=min(np.min(a=lat_lon),np.nanmin(a=lat_dep),np.nanmin(a=dep_lon))
    vmax=max(np.max(a=lat_lon),np.nanmax(a=lat_dep),np.nanmax(a=dep_lon))
    levels=np.linspace(start=vmin,stop=vmax,num=20)

    fig=plt.figure(figsize=(25.6,14.4))
    fig.set_layout_engine(layout="constrained")
    axes=fig.subplots(nrows=2,ncols=2,width_ratios=[2,1],height_ratios=[2,1])
    fig.delaxes(ax=axes[1,1])

    ax=axes[0,0]
    quadcontourset=ax.contourf(lon,lat,lat_lon,levels=levels,cmap="rainbow")
    ax.contour(lon_topo,lat_topo,topo,colors="black",linestyles="solid")
    ax.axvline(x=lon[idx_loc[0]],linestyle="dashed",color="black")
    ax.axhline(y=lat[idx_loc[1]],linestyle="dashed",color="black")
    ax.scatter(x=df["stlo"],y=df["stla"],color="black",marker="^")
    ax.scatter(x=crater[0],y=crater[1],color="black",marker="*")

    ax.set_xlim(left=lon_plt_min,right=lon_plt_max)
    ax.set_ylim(bottom=lat_plt_min,top=lat_plt_max)
    ax.tick_params(labelbottom=False,labelleft=False)
    ax.grid(visible=True)

    ax=axes[0,1]
    ax.contourf(dep,lat,lat_dep,levels=levels,cmap="rainbow")
    ax.axhline(y=lat[idx_loc[1]],linestyle="dashed",color="black")
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
    ax.contourf(lon,dep,dep_lon,levels=levels,cmap="rainbow")
    ax.axvline(x=lon[idx_loc[0]],linestyle="dashed",color="black")
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

    fig.suptitle(t=sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+sys.argv[5])
    fig.savefig(fname=path_save)

def main():
    global lat_lon,lat_dep,dep_lon

    beamform()

    init_topo()

    lat_lon=beam[np.arange(stop=num[0])[:,None],np.arange(stop=num[1]),idx_topo].T
    lat_dep=beam[idx_loc[0]]
    dep_lon=beam[:,idx_loc[1]].T

    plot_3d()

main()
