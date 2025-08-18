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
half_sub=len_sub//2
npts_sub=len_sub*sampling_rate
npts_half_sub=npts_sub//2
lag=np.arange(start=-npts_half_sub,stop=npts_half_sub)/sampling_rate

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
num=vel["num"] #(lon,lat,dep)

#travel time difference
idx_dtt=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/idx_dtt.npz")
idx_dtt=idx_dtt["idx_dtt"] #(ntriu,nnode)

#cross correlation
xcorr=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/xcorr/"+sys.argv[2]+"."+mdl_t.name+".npz")
xcorr=xcorr["xcorr"][idx_avg] #(ntriu,npts_sub)

#location
loc=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/loc/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+mdl_t.name+".npz")
idx_loc=loc["idx_loc"][:,idx_avg] #(3,)

#station
df=pd.read_csv(filepath_or_buffer=mdl_x.info_sta,sep=" ",names=["stnm","stla","stlo","stel"])
df["stdp"]=-0.001*df["stel"]
nsta=len(df)
stnm_pairs=[]
for i in range(nsta):
    for j in range(i+1,nsta):
        stnm_pairs.append(df.loc[i,"stnm"]+"-"+df.loc[j,"stnm"])
ntriu=nsta*(nsta-1)//2

#station pair choice
choice=np.ones(shape=ntriu,dtype=bool)
for i in mdl_b.exclude:
    idx_pair=stnm_pairs.index(i)
    choice[idx_pair]=False

#bandpass filter
if bool_fltr:
    sos_butter=butter(N=4,Wn=[3,6],btype="band",fs=sampling_rate,output="sos")

#array
beam=np.empty(shape=npts_sub)

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/figs/xcf"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+dt.strftime(format="%Y-%m-%dT%H%M%S")+".png"
if isfile(path=path_save):
    print("figure already exists")
    exit()

nstack=0
for i in range(ntriu):
    if np.isnan(xcorr[i,0]) or not choice[i]:
        continue

    if bool_fltr:
        xcorr[i,:]=sosfilt(sos=sos_butter,x=xcorr[i])
        xcorr[i,:]=sosfilt(sos=sos_butter,x=xcorr[i,::-1])[::-1]
    xcorr[i,:]=abs(hilbert(x=xcorr[i]))
    xcorr[i,:]=gaussian_filter1d(input=xcorr[i],sigma=mdl_b.sigma)
    if mdl_b.normalize:
        xcorr[i,:]/=max(xcorr[i])
    else:
        if not nstack:
            norm=max(xcorr[i])
        xcorr[i,:]/=norm

    nstack+=1

def plot_xcf(idx_loc,ax):
    global beam

    idx_node=np.ravel_multi_index(multi_index=idx_loc,dims=num)
    beam[:]=0
    ax.axvline(x=0,color="red")

    j=0
    for i in range(ntriu):
        if np.isnan(xcorr[i,0]) or not choice[i]:
            continue

        idx_shift_dtt=idx_dtt[i,idx_node]-npts_half_sub
        if 0<idx_shift_dtt:
            beam[:-idx_shift_dtt]+=xcorr[i,idx_shift_dtt:]
        elif idx_shift_dtt<0:
            beam[-idx_shift_dtt:]+=xcorr[i,:idx_shift_dtt]
        else: #idx_shift_dtt==0
            beam+=xcorr[i]

        ax.plot(lag-idx_shift_dtt/sampling_rate,j+0.75*xcorr[i],color="black")
        ax.annotate(text=stnm_pairs[i],xy=(-half_sub,j),xytext=(0,5),textcoords="offset points")

        j+=1

    beam/=nstack
    ax.plot(lag,j+beam,color="red")

    ax.grid(visible=True,axis="x")
    ax.set_xlim(left=-half_sub,right=half_sub)
    ax.set_xlabel(xlabel="lag time (s)")
    ax.set_yticks(ticks=[])
    ax.set_ylim(bottom=-1,top=nstack+2)

fig=plt.figure(figsize=(10,14.4))
fig.set_layout_engine(layout="constrained")
axes=fig.subplots(ncols=2)

plot_xcf(idx_loc=idx_loc,ax=axes[0])

fig.suptitle(t=sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+sys.argv[5])
fig.savefig(fname=path_save)
