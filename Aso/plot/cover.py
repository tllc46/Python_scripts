from datetime import datetime

import numpy as np
from matplotlib.ticker import FixedLocator,FixedFormatter
import matplotlib.pyplot as plt

import pandas as pd

df=pd.read_csv(filepath_or_buffer="station_metadata",sep=" ",names=["netwk","stnm","stla","stlo","status"])
df=df[df["status"]=="O"]
nsta=len(df)

#fig=plt.figure(figsize=(8.27,4.65),layout="constrained")
fig=plt.figure(figsize=(25.6,14.4),layout="constrained")
ax=fig.add_subplot()

stnms_formatter=[]
for i in range(nsta):
    netwk=df.iloc[i,0]
    stnm=df.iloc[i,1]
    df_cover=pd.read_csv(filepath_or_buffer="/home/tllc46/48NAS1/tllc46/Aso/avail/entire/"+netwk+"."+stnm+".U.csv",sep=" ",names=["start","end"],parse_dates=[0,1])
    df_cover["width"]=df_cover["end"]-df_cover["start"]
    stnms_formatter.append(stnm)
    if netwk=="N":
        facecolors="green"
    elif netwk=="V":
        facecolors="red"
    ax.broken_barh(xranges=df_cover[["start","width"]].to_numpy(),yrange=(i+1-0.25,0.5),facecolors=facecolors)

ax.set_xlim(left=np.datetime64("2014-04-01"),right=np.datetime64("2015-06-01"))
ax.set_xlabel(xlabel="JST date")
ax.yaxis.set_major_locator(locator=FixedLocator(locs=range(1,nsta+1)))
ax.yaxis.set_major_formatter(formatter=FixedFormatter(seq=stnms_formatter))
ax.grid(visible=True,axis="x")
ax.invert_yaxis()

fig.savefig(fname="/home/tllc46/48NAS1/tllc46/temp_fig/foo106.png")
