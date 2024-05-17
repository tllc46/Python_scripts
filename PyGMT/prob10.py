import pandas as pd
import pygmt

      #event longitude,event latitude,stike,dip,rake,magnitude,plot longitude,plot latitude,text
data=[[129.19,35.77,15,27,71,-179,5.1,128.4,35.5,"2016-09-12T19:44:32"],
      [129.19,35.76,15,24,70,171,5.8,128.5,36.1,"2016-09-12T20:32:54"],
      [129.37,36.11,7,214,51,128,5.4,129.1,36.7,"2017-11-15T14:29:31"]]
focal_mechanism=pd.DataFrame(data=data,columns=["longitude","latitude","depth","strike","dip","rake","magnitude","plot_longitude","plot_latitude","event_name"])

fig=pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="D")
fig.coast(frame=["WeSn","0.5"],land="gray",projection="merc/5c",region=[128,130,35,37])
fig.plot(data="faults.txt",pen="darkolivegreen")
fig.meca(spec=focal_mechanism,scale="0.75c",offset="+gred+s0.1c",compressionfill="blue")
fig.savefig(fname="prob10.png")
