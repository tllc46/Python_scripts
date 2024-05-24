import math
import pygmt

elevation=20

depth_max=60 #km
depth_scale=0.4 #cm/km

y_shift=depth_max*depth_scale*math.cos(math.radians(elevation)) #cm

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D",MAP_FRAME_TYPE="plain")
#depth plot
grid_15=pygmt.xyz2grd(data="depth_15km.txt",spacing=["469+n","325+n"],region=[125.62,127.5115,32.8565,33.95],convention="LB")
grid_30=pygmt.xyz2grd(data="depth_30km.txt",spacing=["469+n","325+n"],region=[125.62,127.5115,32.8565,33.95],convention="LB")
grid_45=pygmt.xyz2grd(data="depth_45km.txt",spacing=["469+n","325+n"],region=[125.62,127.5115,32.8565,33.95],convention="LB")
grid_60=pygmt.xyz2grd(data="depth_60km.txt",spacing=["469+n","325+n"],region=[125.62,127.5115,32.8565,33.95],convention="LB")
pygmt.makecpt(cmap="polar",reverse=True,series=[-0.3,0.3])
fig.grdimage(grid=grid_60,frame="lrbt",projection=["merc/10c",f"z/-{depth_scale}c"],perspective=[150,elevation])
fig.coast(projection="z",shorelines=True,perspective=True)
fig.grdimage(grid=grid_45,frame="lrbt",projection="z",perspective=[150,elevation,-15])
fig.coast(projection="z",shorelines=True,perspective=True)
fig.grdimage(grid=grid_30,frame="lrbt",projection="z",perspective=[150,elevation,-30])
fig.coast(projection="z",shorelines=True,perspective=True)
fig.grdimage(grid=grid_15,frame="lrbt",projection="z",perspective=[150,elevation,-45])
fig.coast(projection="z",shorelines=True,perspective=True)
fig.basemap(frame=["lrbtZ+b","z15","z+ldepth (km)"],region=[125.62,127.5115,32.8565,33.95,0,depth_max],perspective=[150,elevation])
fig.colorbar(frame=[0.1,"x+ldVp (km/s)"],position="JCB+w8c")

#elevation plot
grid=pygmt.datasets.load_earth_relief(resolution="01m",region=[125.6,127.6,32.8,33.95])
fig.shift_origin(yshift=f"{y_shift}c")
fig.grdview(grid=0.001*grid,cmap="geo",shading="+d",projection="merc/10c",zscale="1c",plane=0,surftype="s",region=[125.62,127.5115,32.8565,33.95,0,2],perspective=True) #must provide projection explicitly. colorbar initializes projection by -Jx1i
fig.basemap(frame=["WrStZ","a0.5f0.1","z1","z+lelevation (km)"],perspective=True)
fig.savefig(fname="prob8.png")
