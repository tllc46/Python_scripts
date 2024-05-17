import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="KMA_events.txt",sep=" ")
df=df[["longitude","latitude","depth","magnitude"]]
df["magnitude"]*=0.05

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
region=pygmt.info(data=df[["longitude","latitude"]],spacing=1)
color_range=pygmt.info(data=df["depth"],nearest_multiple=5)[:2]
pygmt.makecpt(cmap="inferno",reverse=True,series=color_range)
fig.coast(frame=["WeSn",2],land="gray",projection="merc/1.5c",region=region)
fig.plot(data="faults3.txt",pen="darkolivegreen")
fig.plot(data=df,cmap=True,style="cc",pen="+cl")
fig.plot(x=[129.05,129.9],y=[35.55,37],pen="blue")
fig.text(x=129.9,y=37.1,text="prob11",justify="LB")
fig.colorbar(frame=[5,"x+ldepth(km)"],position="JCB+w5c")
fig.savefig(fname="prob9.png")
