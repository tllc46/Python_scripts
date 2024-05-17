import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="stations.txt",sep=" ")

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D",MAP_FRAME_TYPE="plain")
grid=pygmt.xyz2grd(data="depth_30km.txt",spacing=["469+n","325+n"],region=[125.62,127.5115,32.8565,33.95],convention="LB")
pygmt.makecpt(cmap="polar",reverse=True,series=[-0.2,0.2])
fig.grdimage(grid=grid,frame=["WeSn","a0.5f0.1"],projection="merc/10c")
fig.coast(shorelines="1p")
fig.plot(data=df[["longitude","latitude"]],fill="white",style="s0.15c",pen="0.75p")
fig.plot(x=[125.62,127.5115],y=[33.35,33.35],pen=",,-")
fig.text(position="RM",text="prob7",justify="RT",offset="j0.5c/0.7c")
fig.colorbar(frame=[0.1,"x+ldVp (km/s)"],position="JCB+w8c")
fig.savefig(fname="prob6.png")
