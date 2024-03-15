from obspy import read
from obspy.core import AttribDict
from obspy.signal.array_analysis import get_geometry,get_timeshift
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

import numpy as np
import matplotlib.pyplot as plt

minslow=10    #minimum slowness
maxslow=20   #maximum slowness
nslow=100     #no. of slowness intervals
tbegin=20   #begin time of timewindow
tend=50     #end time of timewindow

def vesp(slow):
    #convert slowness unit s/° to s/km
    sx=-slow*np.sin(baz*np.pi/180)*180/(np.pi*6371)
    sy=-slow*np.cos(baz*np.pi/180)*180/(np.pi*6371)

    #calculate time shifts
    timeshift=get_timeshift(geometry=geometry,sll_x=sx,sll_y=sy,sl_s=0,grdpts_x=1,grdpts_y=1).squeeze()
    timeshift=delta*np.rint(timeshift/delta)

    st_stack=st.copy()

    #moveout
    for i in range(len(st_stack)):
        st_stack[i].stats.starttime-=timeshift[i]
    
    #common time window
    st_stack.trim(starttime=st[0].stats.starttime+tbegin,endtime=st[0].stats.starttime+tend,pad=True,fill_value=0)

    #stacking
    st_stack.stack(stack_type=("root",4))

    return st_stack[0].data

npvesp=np.frompyfunc(vesp,1,1)

#main
#all trace must align by event time(iztype=io)
st=read(pathname_or_url="filtdata/24/Z/*",format="SAC",byteorder="little")

for i in range(len(st)):
    st[i].stats.coordinates=AttribDict({"latitude":st[i].stats.sac.stla,"longitude":st[i].stats.sac.stlo,"elevation":-12345})

#array geometry
geometry=get_geometry(stream=st,return_center=True)
center=geometry[-1]
geometry=geometry[:-1]

#back azimuth
_,_,baz=gps2dist_azimuth(st[0].stats.sac.evla,st[0].stats.sac.evlo,center[1],center[0])

delta=st[0].stats.delta
time=np.arange(tbegin,tend+delta,delta)
slow=np.linspace(minslow,maxslow,nslow+1)

#beamforming
beam=npvesp(slow)
beam=np.stack(beam)
beam/=beam.max()
time,slow=np.meshgrid(time,slow)

#theorital travel time, ray parameter
model=TauPyModel(model="ak135")
arrivals=model.get_travel_times_geo(st[0].stats.sac.evdp,st[0].stats.sac.evla,st[0].stats.sac.evlo,center[1],center[0],["P","Pn"])

plt.set_cmap("seismic")
plt.pcolormesh(time,slow,beam,vmin=-1,vmax=1)
plt.scatter(arrivals[1].time,arrivals[1].ray_param_sec_degree,c="black",s=240,marker="+")
plt.annotate(arrivals[1].name,(arrivals[1].time,arrivals[1].ray_param_sec_degree),color="black",xytext=(5,5),textcoords="offset points")
plt.scatter(arrivals[2].time,arrivals[2].ray_param_sec_degree,c="black",s=240,marker="+")
plt.annotate(arrivals[2].name,(arrivals[2].time,arrivals[2].ray_param_sec_degree),color="black",xytext=(5,5),textcoords="offset points")
plt.xlabel("time(s)")
plt.ylabel("slowness(s/°)")
plt.show()
