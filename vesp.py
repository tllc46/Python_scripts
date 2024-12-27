from obspy import read
from obspy.core import AttribDict
from obspy.signal.array_analysis import get_geometry,get_timeshift
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

import numpy as np
import matplotlib.pyplot as plt

minslow=10 #최소 slowness
maxslow=20 #최대 slowness
nslow=100 #slowness 간격 개수
tbegin=20 #stack 구간 시작 시각
tend=50 #stack 구간 끝 시각

def vesp(slow):
    #slowness 단위를 s/degree->s/km 변환
    sx=-slow*np.sin(np.radians(baz))*180/(np.pi*6371)
    sy=-slow*np.cos(np.radians(baz))*180/(np.pi*6371)

    #time shift 계산
    timeshift=get_timeshift(geometry=geometry,sll_x=sx,sll_y=sy,sl_s=0,grdpts_x=1,grdpts_y=1).squeeze()
    timeshift=delta*np.rint(timeshift/delta)

    st_stack=st.copy()

    #moveout
    for i in range(len(st_stack)):
        st_stack[i].stats.starttime-=timeshift[i]
    
    #stack할 구간으로 자르기
    st_stack.trim(starttime=st[0].stats.starttime+tbegin,endtime=st[0].stats.starttime+tend,pad=True,fill_value=0)

    #stacking
    st_stack.stack(stack_type=("root",4))

    return st_stack[0].data

npvesp=np.frompyfunc(vesp,nin=1,nout=1)

#main
#모든 trace는 event 발생 시각으로 정렬되어 있어야 한다(iztype=io)
st=read(pathname_or_url="filtdata/24/Z/*",format="SAC",byteorder="little")

for i in range(len(st)):
    st[i].stats.coordinates=AttribDict({"latitude":st[i].stats.sac.stla,"longitude":st[i].stats.sac.stlo,"elevation":-12345})

#지진계 중심으로부터 거리 vector
geometry=get_geometry(stream=st,return_center=True)
center=geometry[-1]
geometry=geometry[:-1]

#back azimuth
_,_,baz=gps2dist_azimuth(lat1=st[0].stats.sac.evla,lon1=st[0].stats.sac.evlo,lat2=center[1],lon2=center[0])

delta=st[0].stats.delta
time=np.arange(stop=int((tend-tbegin)/delta)+1)*delta+tbegin
slow=np.arange(stop=nslow+1)*(maxslow-minslow)/nslow+minslow

#beamforming
beam=npvesp(slow)
beam=np.stack(arrays=beam)
beam/=beam.max()

#theorital travel time, ray parameter
model=TauPyModel(model="ak135")
arrivals=model.get_travel_times_geo(source_depth_in_km=st[0].stats.sac.evdp,source_latitude_in_deg=st[0].stats.sac.evla,source_longitude_in_deg=st[0].stats.sac.evlo,receiver_latitude_in_deg=center[1],receiver_longitude_in_deg=center[0],phase_list=["P","Pn"])

fig=plt.figure()
ax=fig.subplots()
ax.pcolormesh(time,slow,beam,cmap="seismic",vmin=-1,vmax=1)
ax.scatter(x=arrivals[1].time,y=arrivals[1].ray_param_sec_degree,s=240,c="black",marker="+")
ax.annotate(text=arrivals[1].name,xy=(arrivals[1].time,arrivals[1].ray_param_sec_degree),xytext=(5,5),textcoords="offset points",color="black")
ax.scatter(x=arrivals[2].time,y=arrivals[2].ray_param_sec_degree,s=240,c="black",marker="+")
ax.annotate(text=arrivals[2].name,xy=(arrivals[2].time,arrivals[2].ray_param_sec_degree),xytext=(5,5),textcoords="offset points",color="black")
ax.set_xlabel(xlabel="time(s)")
ax.set_ylabel(ylabel="slowness(s/°)")
plt.show()
