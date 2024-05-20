import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="stations.txt",sep=" ")

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D",MAP_FRAME_TYPE="plain")
grid=pygmt.xyz2grd(data="depth_30km.txt",spacing=["469+n","325+n"],region=[125.62,127.5115,32.8565,33.95],convention="LB")
series=pygmt.grdinfo(grid=grid,nearest_multiple="0.05+s")
series=series[2:-1].split(sep="/")[:2]
pygmt.makecpt(cmap="polar",reverse=True,series=series)
fig.grdimage(grid=grid,frame=["WeSn","a0.5f0.1"],projection="merc/5c")
fig.coast(shorelines=True)
fig.plot(data=df[["longitude","latitude"]],fill="white",style="s0.12c",pen=True)
fig.plot(x=[125.62,127.5115],y=[33.35,33.35],pen="dashed")
fig.text(x=127.5115,y=33.35,text="prob7",justify="RT",offset="j0.1c")
fig.colorbar(frame=[0.1,"x+ldVp (km/s)"],position="JCB+w4c")
fig.savefig(fname="prob6.png")
