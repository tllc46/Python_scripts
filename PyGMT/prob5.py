import numpy as np
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
import pygmt

#calculate baz
df=pd.read_csv(filepath_or_buffer="traveltime_residuals.txt",sep=" ")
longitude_center,latitude_center=np.average(a=df[["longitude","latitude"]],axis=0)
gca,az,baz=gps2dist_azimuth(lat1=-12.0989,lon1=166.5894,lat2=latitude_center,lon2=longitude_center)
df_O=df.loc[df["status"].isin(["O"]),["longitude","latitude","residual"]]
df_O["size"]=2*abs(df_O.loc[:,"residual"])+0.1

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D",MAP_FRAME_TYPE="plain")
#main plot
fig.coast(frame=["WeSn","xa0.2-0.1f0.1","y0.15"],projection="merc/10c",region=[126.1,127,33.15,33.6],shorelines=True)
color_range=pygmt.info(data=df["residual"],nearest_multiple=0.05)[:2]
pygmt.makecpt(cmap="polar",series=color_range)
fig.plot(data=df_O,cmap=True,style="cc",pen=True)
fig.plot(data=df.loc[df["status"].isin(["X"]),["longitude","latitude"]],style="+0.15c")
fig.text(x=df["longitude"],y=df["latitude"],text=df["name"],font="5p",justify="RB",offset="j0.1c",fill="white@40")
fig.colorbar(frame=[0.05,"x+lresidual (s)"],position="JCB+w4c")

#inset plot
with fig.inset(position="jRB+o0.1c",projection="xy/1.5c",region=[0,1,0,1]):
    fig.plot(data=[[0.5,0.5,baz,0.7]],style="V0.2c+e")
    fig.text(position="CT",text="Baz.",font="7p",justify="CT",offset="0.2c/-0.2c")
fig.savefig(fname="prob5.png")
