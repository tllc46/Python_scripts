from io import StringIO
import pandas as pd
import pygmt

df=pd.read_csv(filepath_or_buffer="KMA_events.txt",sep=" ")

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
region=pygmt.info(data=df[["longitude","latitude"]],spacing=1)
series=pygmt.info(data=df["depth"],nearest_multiple=5)[:2]
pygmt.makecpt(cmap="inferno",reverse=True,series=series)
fig.coast(frame=["WeSn",2],land="gray",projection="merc/1.5c",region=region)
fig.plot(data="faults.txt",pen="darkolivegreen")
fig.plot(x=df["longitude"],y=df["latitude"],size=0.05*df["magnitude"],cmap=True,fill=df["depth"],style="cc",pen="+cl")
fig.plot(x=[129.05,129.9],y=[35.55,37],pen="blue")
fig.text(x=129.9,y=37.1,text="prob11",justify="LB")
fig.colorbar(frame=[5,"x+ldepth(km)"],position="JCB+w5c")

#legend
legend_spec=StringIO(initial_value="""
H 7p magnitude
S 0.4c c 0.1c - ,, 0.7c 2
G 0.05c
S 0.4c c 0.15c - ,, 0.7c 3
G 0.05c
S 0.4c c 0.2c - ,, 0.7c 4
G 0.05c
S 0.4c c 0.25c - ,, 0.7c 5
G 0.05c
S 0.4c c 0.3c - ,, 0.7c 6
""")
with pygmt.config(FONT_ANNOT_PRIMARY="7p"):
    fig.legend(spec=legend_spec,position="JRM+w1.5c+o0.2c/0c")

fig.savefig(fname="prob9.png")
