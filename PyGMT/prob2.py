import numpy as np
import pygmt

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
#main plot
grid=pygmt.datasets.load_earth_relief(resolution="03s",region=[126,127,33.15,33.6])
fig.grdimage(grid=grid,frame=["WeSn",0.2],cmap="geo",shading="+a120",projection="merc/10c")
fig.colorbar(frame=[500,"x+lelevation","y+lm"],position="JCB+w4c")

#inset plot
with fig.inset(position="jLT+o0.1c",box="+gwhite+p0.4p",projection="merc/0.28c",region=[125,132,32,38]):
    fig.coast(land="gray")
    fig.plot(data=np.array([[126,33.15,127,33.6]]),style="r+s",pen="blue")
with fig.inset(position="jLT+o0.1c",box="+p0.4p",projection="merc/0.28c",region=[125,132,32,38]):
    pass
fig.savefig(fname="prob2.png")
