#usage
#python main_bp.py t01 x01 v01 b01 2014

import sys
from os import makedirs
from os.path import isdir,isfile
from importlib import import_module

import numpy as np
from scipy.signal import butter,sosfilt,hilbert
from scipy.ndimage import gaussian_filter1d
import pandas as pd

bool_fltr=False

sys.path.append("/home/tllc46/Aso/xcorr")

mdl_x=import_module(name=sys.argv[2])
mdl_b=import_module(name=sys.argv[4])

#sampling rate
sampling_rate=100 #[Hz]

#sub window
len_sub=48 #[s]
npts_sub=len_sub*sampling_rate

#velocity
vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/vel.npz")
idx_flat=vel["idx"].flatten() #(lon,lat,dep) → (nnode,)
num=vel["num"] #(lon,lat,dep)
nnode=np.prod(a=num)

#travel time difference
idx_dtt=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[3]+"/idx_dtt.npz")
idx_dtt=idx_dtt["idx_dtt"] #(ntriu,nnode)

#cross correlation
xcorr=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/xcorr/"+sys.argv[2]+"."+sys.argv[5]+".npz")
xcorr=xcorr["xcorr"] #(navg,ntriu,npts_sub)
navg=xcorr.shape[0]

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
xcorr_pair=np.empty(shape=npts_sub)
beam=np.empty(shape=nnode)
idx_loc=np.empty(shape=(3,navg),dtype=int)
beam_max=np.empty(shape=navg)

#saving directory
path_save="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/loc"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[2]+"."+sys.argv[3]+"."+sys.argv[4]+"."+sys.argv[5]+".npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

def beamform():
    global xcorr_pair,beam
    nstack=0
    beam[:]=0

    for i in range(ntriu):
        xcorr_pair[:]=xcorr[idx_avg,i]
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
    beam_max[idx_avg]=beam[idx_max]
    idx_loc[:,idx_avg]=np.unravel_index(indices=idx_max,shape=num)

def main():
    global idx_avg

    for idx_avg in range(navg):
        beamform()

    np.savez(file=path_save,idx_loc=idx_loc,beam_max=beam_max)

main()
