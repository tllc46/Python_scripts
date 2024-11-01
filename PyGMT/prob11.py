import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="KMA_events.txt",sep=" ")
df_big=df.loc[df["magnitude"].isin([5.1,5.8,5.4])]

depth_scale=0.1 #cm/km
depth_max=25 #km
y_shift=depth_max*depth_scale #cm

fig=pygmt.Figure()
pygmt.config(FONT_ANNOT_PRIMARY="6p",FONT_LABEL="8p",FORMAT_GEO_MAP="D")
#map plot
fig.coast(frame="ag",land="gray",projection="omerc/129.05/35.55/129.9/37/6c",region=[-0.05,1.65,-0.2,0.2])
fig.plot(x=[129.05,129.9],y=[35.55,37],pen="0.75p,blue")
fig.plot(x=df["longitude"],y=df["latitude"],size=0.05*df["magnitude"],style="cc",pen="red")
fig.plot(x=df_big["longitude"],y=df_big["latitude"],size=0.05*df_big["magnitude"],style="cc",pen="black")

#height plot
grid=pygmt.datasets.load_earth_relief(resolution="01m",region=[129,130,35,38])
track=pygmt.grdtrack(grid=grid,profile="129.05/35.55/129.9/37+d+i0.01d")
track=track[[2,3]]
track[3]*=0.001
region=pygmt.info(data=track,spacing=[0,0.2])
fig.shift_origin(xshift="0.3c",yshift="-2c")
fig.plot(x=region[:2],y=[0,0],frame=["lEt","y1","y+lheight(km)"],fill="lightblue",projection="xy/6c/0.5c",close="+yb",region=region)
fig.plot(x=track[2],y=track[3],fill="black",close="+yb")

#depth plot
profile=pygmt.project(data=df[["longitude","latitude","depth","magnitude"]],center=[129.05,35.55],endpoint=[129.9,37],convention="pz",length="w",width=[-0.2,0.2])
profile_big=profile.loc[profile[2].isin([5.1,5.8,5.4])]
region[2:]=pygmt.info(data=profile[1],nearest_multiple=5)[:2]
fig.shift_origin(yshift=f"-{y_shift}c")
fig.plot(x=profile[0],y=profile[1],size=0.05*profile[2],frame=["WrS","x0.5","x+ldistance(@.)","ya10f5","y+ldepth(km)"],projection=f"xy/6c/-{depth_scale}c",region=region,style="cc",pen="red")
fig.plot(x=profile_big[0],y=profile_big[1],size=0.05*profile_big[2],style="cc",pen="black")
fig.savefig(fname="prob11.png")
