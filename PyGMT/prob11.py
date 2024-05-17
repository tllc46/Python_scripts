import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="KMA_events.txt",sep=" ")

lon_1=129.05
lat_1=35.55
lon_2=129.9
lat_2=37
lon_center=0.5*(lon_1+lon_2)
lat_center=0.5*(lat_1+lat_2)

depth_scale=0.1 #cm/km
depth_max=25 #km
y_shift=depth_max*depth_scale #cm

fig=pygmt.Figure()
pygmt.config(FONT_ANNOT_PRIMARY="6p",FONT_LABEL="8p",FORMAT_GEO_MAP="D")
#map plot
fig.coast(frame="ag",land="gray",projection=f"omerc/{lon_center}/{lat_center}/{lon_2}/{lat_2}/6c",region=[-0.85,0.85,-0.2,0.2])
fig.plot(x=[lon_1,lon_2],y=[lat_1,lat_2],pen="0.75p,blue")
fig.plot(x=df["longitude"],y=df["latitude"],size=0.05*df["magnitude"],style="cc",pen="red")
fig.plot(x=df.loc[df["magnitude"].isin([5.1,5.8,5.4]),"longitude"],y=df.loc[df["magnitude"].isin([5.1,5.8,5.4]),"latitude"],size=0.05*df.loc[df["magnitude"].isin([5.1,5.8,5.4]),"magnitude"],style="cc",pen="black")

#elevation plot
grid=pygmt.datasets.load_earth_relief(resolution="01m",region=[129,130,35,38])
track=pygmt.grdtrack(grid=grid,profile=f"{lon_1}/{lat_1}/{lon_2}/{lat_2}+d+i0.01d")
height_range=pygmt.info(data=0.001*track[3],nearest_multiple=0.2)
fig.shift_origin(xshift="0.3c",yshift="-2c")
fig.plot(x=[0,track.iloc[-1,2]],y=[0,0],frame=["lEt","y1","y+lheight(km)"],fill="lightblue",projection="xy/6c/0.5c",close="+yb",region=[0,track.iloc[-1,2],height_range[0],height_range[1]])
fig.plot(x=track[2],y=0.001*track[3],fill="black",close="+yb")

#depth plot
profile=pygmt.project(data=df.loc[:,["longitude","latitude","depth","magnitude"]],center=[lon_1,lat_1],endpoint=[lon_2,lat_2],width=[-0.2,0.2])
fig.shift_origin(yshift=f"-{y_shift}c")
fig.plot(x=profile[4],y=profile[2],size=0.05*profile[3],frame=["WrS","x0.5","x+ldistance(@.)","ya10f5","y+ldepth(km)"],projection=f"xy/6c/-{depth_scale}c",region=[0,track.iloc[-1,2],0,depth_max],style="cc",pen="red")
fig.plot(x=profile.loc[profile[3].isin([5.1,5.8,5.4]),4],y=profile.loc[profile[3].isin([5.1,5.8,5.4]),2],size=0.05*profile.loc[profile[3].isin([5.1,5.8,5.4]),3],style="cc",pen="black")
fig.savefig(fname="prob11.png")
