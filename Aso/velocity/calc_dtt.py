#usage
#v06_main.py → set_grid.py → calc_tt.py → calc_dtt.py
#python calc_dtt.py v06

import sys
from os.path import isfile

import numpy as np
import pandas as pd

sampling_rate=100 #[Hz]

path_save="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/idx_dtt.npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

df=pd.read_csv(filepath_or_buffer="center",sep=" ",names=["stnm","stla","stlo","stel"])
nsta=len(df)
ntriu=nsta*(nsta-1)//2

grid=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/grid.npz")
num=grid["num"] #(lon,lat,dep)
nnode=np.prod(a=num)

travel_times=np.empty(shape=(nsta,nnode))
diff_travel_times=np.empty(shape=nnode)
idx_dtt=np.empty(shape=(ntriu,nnode),dtype=int)

for i in range(nsta):
    stnm=df.loc[i,"stnm"]
    travel_times_flat=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/vel/"+sys.argv[1]+"/"+stnm+".npz") #(lon,lat,dep)
    travel_times[i,:]=travel_times_flat["travel_times"].flatten() #(nnode,)

for i in range(nsta):
    for j in range(i+1,nsta):
        idx_triu=i*(nsta-1)-i*(i+1)//2+(j-1)
        diff_travel_times[:]=travel_times[i]-travel_times[j]
        idx_dtt[idx_triu,:]=np.round(a=sampling_rate*diff_travel_times).astype(dtype=int)

np.savez(file=path_save,idx_dtt=idx_dtt)
