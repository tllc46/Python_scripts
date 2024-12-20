import math
import pygmt

x_min=125.62 #degree
x_max=127.5115 #degree
x_center=0.5*(x_min+x_max) #degree

depth_min=-65.4 #km
depth_scale=0.075 #cm/km
depth_radius_val=6371.007 #km
radius_len=depth_radius_val*depth_scale #cm

height_min=-0.5 #km
height_max=2 #km
height_scale=0.5 #cm/km
height_radius_val=radius_len/height_scale #km

y_shift=(height_min*height_scale-depth_min*depth_scale)*math.cos(math.radians(x_max-x_center)) #cm

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
#depth plot
grid=pygmt.xyz2grd(data="latitude_33.35.txt",spacing=["703+n","145+n"],region=[x_min,x_max,depth_min+depth_radius_val,0.6+depth_radius_val],convention="LB")
series=pygmt.grdinfo(grid=grid,nearest_multiple="0.05+s")
series=series[2:-1].split(sep="/")[:2]
pygmt.makecpt(cmap="polar",reverse=True,series=series)
fig.grdimage(grid=grid,frame=["WrS","xa1f0.2","y5","y+lDepth(km)"],projection=f"polar/{depth_scale}c+a+t{x_center}+z")
fig.text(position="CB",text="Longitude",font="12p",justify="CT",offset="j0c/0.5c",no_clip=True)
fig.colorbar(frame=[0.1,"x+ldVp (km/s)"],position="JCB+w6c+o0c/1.5c")

#height plot
grid_relief=pygmt.datasets.load_earth_relief(resolution="01m",region=[125.6,127.6,32.8,33.95])
track=pygmt.grdtrack(grid=grid_relief,profile=f"{x_min}/33.35/{x_max}/33.35+d+i0.01d")
fig.shift_origin(yshift=f"{y_shift}c")
fig.plot(x=[x_min,x_max],y=[height_radius_val,height_radius_val],frame="+gwhite",fill="lightblue",projection=f"polar/{height_scale}c+a+t{x_center}",close="+yb",region=[x_min,x_max,height_radius_val+height_min,height_radius_val+height_max])
fig.plot(x=track[0],y=0.001*track[3]+height_radius_val,fill="black",close="+yb")
fig.basemap(frame=["lEt","y1","y+lHeight(km)"],projection=f"polar/{height_scale}c+a+r{radius_len}+t{x_center}",region=[x_min,x_max,height_min,height_max])
fig.savefig(fname="prob7.png")
