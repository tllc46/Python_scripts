import pygmt

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
grid=pygmt.datasets.load_earth_relief(resolution="01m",region=[124,132,33,43])
fig.grdimage(grid=grid,frame=["WeSn","x2","y2+1"],cmap="geo",projection="merc/1c")
fig.coast(borders=[1,"0.75p"])
fig.colorbar(frame=[2500,"x+lelevation","y+lm"],position="JCB+w4c")
fig.savefig(fname="prob1.png")
