import numpy as np
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
import pygmt

#calculate baz
df=pd.read_csv(filepath_or_buffer="traveltime_residuals.txt",sep=" ")
gca,az,baz=gps2dist_azimuth(lat1=-12.0989,lon1=166.5894,lat2=np.average(df["latitude"]),lon2=np.average(df["longitude"]))

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D",MAP_FRAME_TYPE="plain")
#main plot
fig.coast(frame=["WeSn","xa0.2-0.1f0.1","y0.15"],projection="merc/10c",region=[126.1,127,33.15,33.6],shorelines=True)
color_range=pygmt.info(data=df["residual"],nearest_multiple=0.05)[:2]
pygmt.makecpt(cmap="polar",series=color_range)
fig.plot(x=df.loc[df["status"]=="O","longitude"],y=df.loc[df["status"]=="O","latitude"],size=2*abs(df.loc[df["status"]=="O","residual"])+0.1,cmap=True,fill=df.loc[df["status"]=="O","residual"],style="cc",pen=True)
fig.plot(x=df.loc[df["status"]=="X","longitude"],y=df.loc[df["status"]=="X","latitude"],style="+0.15c")
fig.text(x=df["longitude"],y=df["latitude"],text=df["name"],font="5p",justify="RB",offset="j0.1c",fill="white@40")
fig.colorbar(frame=[0.05,"x+lresidual (s)"],position="JCB+w4c")

#inset plot
with fig.inset(position="jRB+o0.1c",projection="xy/1.5c",region=[0,1,0,1]):
    fig.plot(data=np.array([[0.5,0.5,baz,0.7]]),style="V0.2c+e")
    fig.text(position="CT",text="Baz.",font="7p",justify="CT",offset="0.2c/-0.2c")
fig.savefig(fname="prob5.png")