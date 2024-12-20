import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="stations.txt",sep=" ")
df_KMA=df.loc[df["network"].isin(["KMA"])]
df_KIGAM=df.loc[df["network"].isin(["KIGAM"])]
df_Temp=df.loc[df["network"].isin(["Temp"])]

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
grid=pygmt.datasets.load_earth_relief(resolution="03s",region=[126.1,127,33.15,33.6])
grid_grad=pygmt.grdgradient(grid=grid,azimuth=120,normalize="t")
pygmt.makecpt(cmap="gray",series=[-1,1])
fig.grdimage(grid=grid_grad,frame=["WeSn",0.2],projection="merc/10c")
fig.coast(water="white",shorelines="0.4p")
fig.plot(x=df_KMA["longitude"],y=df_KMA["latitude"],fill="red",style="t0.2c",pen=True,label="KMA")
fig.plot(x=df_KIGAM["longitude"],y=df_KIGAM["latitude"],fill="green",style="t0.2c",pen=True,label="KIGAM+G0.1c")
fig.plot(x=df_Temp["longitude"],y=df_Temp["latitude"],fill="blue",style="t0.2c",pen=True,label="temporary+G0.1c")
fig.text(x=df["longitude"],y=df["latitude"],text=df["name"],font="5p",justify="CT",offset="j0.1c",fill="white@40")
with pygmt.config(FONT_ANNOT_PRIMARY="5p"):
    fig.legend(position="jRB+o0.1c",box="+gwhite+p0.5p")
fig.savefig(fname="prob3.png")
