#usage
#v01_main.py → calc_tt.py → vel_reduce.py → calc_dtt.py
#python calc_dtt.py v01

import sys
from os.path import isfile

import numpy as np
import pandas as pd

sampling_rate=100 #[Hz]
len_sub=48 #[s]
npts_sub=len_sub*sampling_rate

path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/idx_dtt.npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

df=pd.read_csv(filepath_or_buffer="center",sep=" ",names=["stnm","stla","stlo","stel"])
nsta=len(df)
ntriu=nsta*(nsta-1)//2

vel=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/vel.npz")
num=vel["num"] #(lon,lat,dep)
nnode=np.prod(a=num)

travel_times=np.empty(shape=(nsta,nnode))
diff_travel_times=np.empty(shape=nnode)
idx_dtt=np.empty(shape=(ntriu,nnode),dtype=int)

for i in range(nsta):
    stnm=df.loc[i,"stnm"]
    travel_times_flat=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/"+stnm+".npz")
    travel_times[i,:]=travel_times_flat["travel_times"] #(nnode,)

for i in range(nsta):
    for j in range(i+1,nsta):
        idx_triu=i*(nsta-1)-i*(i+1)//2+(j-1)
        diff_travel_times[:]=travel_times[i]-travel_times[j]
        idx_dtt[idx_triu,:]=np.round(a=sampling_rate*diff_travel_times).astype(dtype=int)

np.savez(file=path_save,idx_dtt=idx_dtt)
