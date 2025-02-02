import pandas as pd
import pygmt

radius_scale=0.05 #cm/r
radius_len=90*radius_scale #cm, aeqd's center of projection corresponds to oblique latitude 90

#calculate center of stations
df_stations=pd.read_csv(filepath_or_buffer="stations.txt",sep=" ")
longitude_center,latitude_center=df_stations[["longitude","latitude"]].mean()
df=pd.read_csv(filepath_or_buffer="events.txt",sep=" ")

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="+D")
fig.coast(land="whitesmoke",projection=f"aeqd/{longitude_center}/{latitude_center}/120/{radius_len}c/0",region="g",water="darkgray")
fig.basemap(frame=["xa30g30","yg30"],projection=f"polar/{radius_scale}c+a",region=[0,360,0,120])
fig.text(x=[60,60,60],y=[30,60,90],text=["30@.","60@.","90@."],justify="LT",offset="j0.11c/0c")
fig.plot(x=df["longitude"],y=df["latitude"],fill="white",projection=f"aeqd/{longitude_center}/{latitude_center}/120/{radius_len}c/0",region="g",style="c0.15c",pen=True)
fig.plot(x=166.5894,y=-12.0989,fill="red",style="c0.18c",pen=True)
fig.text(x=166.5894,y=-12.0989,text="prob5",justify="RT",offset="j0.1c")
fig.savefig(fname="prob4.png")
