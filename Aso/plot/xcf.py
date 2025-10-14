#usage
#python plot_xcf.py t01 x02 v01 b04 2016-10-19T03:00:00

import sys
from os import makedirs
from os.path import isdir,isfile
from importlib import import_module
from datetime import datetime

import numpy as np
from scipy.signal import butter,sosfilt,hilbert
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import pandas as pd

bool_fltr=False

sys.path.append("/home/tllc46/Aso/tstep")
sys.path.append("/home/tllc46/Aso/xcorr/params")
sys.path.append("/home/tllc46/Aso/bp_ext")

tstep=int(sys.argv[1][1:])+1
tstep=f"t{tstep:02}"

mdl_t_x=import_module(name=sys.argv[1])
mdl_t_l=import_module(name=tstep)
mdl_x=import_module(name=sys.argv[2])
mdl_b=import_module(name=sys.argv[4])

#sampling rate
sampling_rate=100 #[Hz]

#average window
dt=datetime.strptime(sys.argv[5],"%Y-%m-%dT%H:%M:%S")
mdl_t_x.parse_dt(dt=dt)
mdl_t_l.parse_dt(dt=dt)

#velocity grid
grid=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/grid.npz")
num=grid["num"] #(lon,lat,dep)

#cross correlation
xcorr=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/xcorr/"+sys.argv[2]+"."+mdl_t_x.name+".npz")
xcorr=xcorr["xcorr"][mdl_t_x.idx_avg] #(ntriu,npts_sub)
npts_sub=xcorr.shape[1]

#travel time difference
idx_dtt=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/idx_dtt.npz")
idx_dtt=idx_dtt["idx_dtt"] #(ntriu,nnode)

#sub window
len_sub=npts_sub//sampling_rate
lag=np.arange(start=-(npts_sub//2),stop=npts_sub//2)/sampling_rate
lag_max=15 #[s]

#location
loc=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+tstep+"/loc/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+mdl_t_l.name+".npz")
idx_loc=loc["idx_loc"][:,mdl_t_l.idx_avg] #(3,)

#station
df=pd.read_csv(filepath_or_buffer=mdl_x.info_sta,sep=" ",names=["stnm","stla","stlo","stel"])
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
xcorr_original=np.copy(a=xcorr)
beam=np.empty(shape=npts_sub)

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/figs/xcf"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+dt.strftime(format="%Y-%m-%dT%H%M%S")+".png"

nstack=0
for i in range(ntriu):
    if np.isnan(xcorr[i,0]) or not choice_pair[i]:
        continue

    if bool_fltr:
        xcorr[i,:]=sosfilt(sos=sos_butter,x=xcorr[i])
        xcorr[i,:]=sosfilt(sos=sos_butter,x=xcorr[i,::-1])[::-1]
    xcorr[i,:]=abs(hilbert(x=xcorr[i]))
    xcorr[i,:]=gaussian_filter1d(input=xcorr[i],sigma=mdl_b.sigma)
    if mdl_b.normalize:
        xcorr[i,:]/=max(xcorr[i])
        xcorr_original[i,:]/=max(abs(xcorr_original[i]))
    else:
        if not nstack:
            norm=max(xcorr[i])
            norm_original=max(abs(xcorr_original[i]))
        xcorr[i,:]/=norm
        xcorr_original[i,:]/=norm_original

    nstack+=1

#figure
scale=0.65

def plot_xcf(idx_loc,ax):
    global beam

    idx_node=np.ravel_multi_index(multi_index=idx_loc,dims=num)
    beam[:]=0

    j=0
    for i in range(ntriu):
        if np.isnan(xcorr[i,0]) or not choice_pair[i]:
            continue

        idx_shift_dtt=idx_dtt[i,idx_node]
        if 0<idx_shift_dtt:
            beam[:-idx_shift_dtt]+=xcorr[i,idx_shift_dtt:]
        elif idx_shift_dtt<0:
            beam[-idx_shift_dtt:]+=xcorr[i,:idx_shift_dtt]
        else: #idx_shift_dtt==0
            beam+=xcorr[i]

        ax.plot(lag-idx_shift_dtt/sampling_rate,j+scale*xcorr_original[i],color="blue")
        ax.plot(lag-idx_shift_dtt/sampling_rate,j+scale*xcorr[i],color="black")
        ax.annotate(text=stnm_pairs[i],xy=(-lag_max,j),xytext=(1,5),textcoords="offset points")

        j+=1

    beam/=nstack
    ax.plot(lag,j+scale*beam,color="red")
    ax.axvline(x=0,color="red")

    ax.grid(visible=True,axis="x")
    ax.set_xlim(left=-lag_max,right=lag_max)
    ax.set_xlabel(xlabel="lag time (s)")
    ax.set_yticks(ticks=[])
    ax.set_ylim(bottom=-1,top=nstack+1)

def main():
    #fig=plt.figure(figsize=(8.27,11.69),layout="constrained")
    fig=plt.figure(figsize=(25.6,14.4),layout="constrained")
    axes=fig.subplots(ncols=2)

    plot_xcf(idx_loc=idx_loc,ax=axes[0])
    plot_xcf(idx_loc=[123,120,23],ax=axes[1])

    fig.savefig(fname=path_save)

main()
