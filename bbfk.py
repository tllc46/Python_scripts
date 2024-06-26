from obspy import read
from obspy.core import AttribDict
from obspy.signal.array_analysis import get_geometry
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

import numpy as np
import matplotlib.pyplot as plt

eps=0
exp=1
wavenumber=0.02
ssq=100

def phaseshift(iptr,ncurrent):
    #dbh/phaseshift.c/phaseshift()
    for i in range(nch):
        [x1,x2,x3,y1,y2,w1,w2,p1,p2,q1,q2]=qstates[i]

        for j in range(ncurrent):
            y=0.94167*(y2-x1)+x3
            output1=0.53239*(w2-y)+y2
            y2=y1
            y1=y
            w2=w1
            w1=output1

            p=0.186540*(p2-st[i].data[iptr+j])+x2
            output2=0.7902015*(q2-p)+p2
            p2=p1
            p1=p
            q2=q1
            q1=output2

            x3=x2
            x2=x1
            x1=st[i].data[iptr+j]

            output[i][j]=output1+1j*output2

        qstates[i]=[x1,x2,x3,y1,y2,w1,w2,p1,p2,q1,q2]

def covmat():
    #fks/xbbfk.c/covmat()
    decrate=2

    for iptr in range(0,nsamples,100): #100 samples window로 끊어서 scm 계산한 뒤 합
        ncurrent=min(nsamples-iptr,100)
        phaseshift(iptr=iptr,ncurrent=ncurrent)
        for i in range(nch):
            for j in range(nch):
                x=output[i][1:ncurrent:decrate] #1/2 desample
                y=output[j][1:ncurrent:decrate]
                c=x*np.conjugate(y)
                scm[j][i]+=sum(c)

def normalize():
    trace=np.real(val=np.trace(scm))
    for i in range(nch):
        v=np.real(scm[i][i])
        scale=np.sqrt(trace/(v*nch))
        for j in range(nch):
            scm[i][j]*=scale
            scm[j][i]*=scale

def regularize():
    trace=np.real(val=np.trace(scm))
    offset=trace*eps/nch
    for i in range(nch):
        scm[i][i]+=offset

def mlm():
    #fks/xbbfk.c/eigenanal()
    w,z=np.linalg.eig(a=scm)
    tmp=np.matmul(z,np.diag(v=1/w))
    sinv=np.matmul(tmp,np.conjugate(np.transpose(z)))
    return sinv

def music():
    #fks/xbbfk.c/eigenanal()
    w,z=np.linalg.eig(scm)
    ind=np.argsort(np.real(w))
    ind=ind[:-1]
    sinv=np.matmul(z[:,ind],np.conjugate(np.transpose(a=z[:,ind])))
    return sinv

def fkevalr():
    #fks/xbbfk.c/fkevalr()
    for i in range(ssq):
        for j in range(ssq):
            tmp=np.matmul(np.exp(-1j*(wv[i]*geometry[:,0]+wv[j]*geometry[:,1])),scm)
            fks[i][j]=np.real(val=np.matmul(tmp,np.exp(1j*(wv[i]*geometry[:,0]+wv[j]*geometry[:,1]))))

#main
st=read(pathname_or_url="sacdata/*",format="SAC",byteorder="little")
nch=len(st)

start=[]
end=[]
for i in range(nch):
    start.append(st[i].stats.starttime)
    end.append(st[i].stats.endtime)
start_com=max(start)
end_com=min(end)
st=st.slice(starttime=start_com,endtime=end_com)
nsamples=st[0].stats.npts

for i in range(nch):
    st[i].stats.coordinates=AttribDict({"latitude":st[i].stats.sac.stla,"longitude":st[i].stats.sac.stlo,"elevation":-12345})

geometry=get_geometry(stream=st,return_center=True)
center=geometry[-1]
geometry=geometry[:-1]
_,_,baz=gps2dist_azimuth(lat1=st[0].stats.sac.evla,lon1=st[0].stats.sac.evlo,lat2=center[1],lon2=center[0])
baz=np.radians(baz)

output=np.empty(shape=(nch,100),dtype=complex)
scm=np.zeros(shape=(nch,nch),dtype=complex)
fks=np.empty(shape=(ssq,ssq))
qstates=np.zeros(shape=(nch,11))
k=np.arange(stop=ssq)*2*wavenumber/(ssq-1)+wavenumber
wv=2*np.pi*k

covmat()
scm=music()
fkevalr()

model=TauPyModel(model="ak135")
arrivals=model.get_travel_times_geo(source_depth_in_km=st[0].stats.sac.evdp,source_latitude_in_deg=st[0].stats.sac.evla,source_longitude_in_deg=st[0].stats.sac.evlo,receiver_latitude_in_deg=center[1],receiver_longitude_in_deg=center[0],phase_list=["P"])
p=arrivals[0].ray_param_sec_degree
p*=0.05/111.19 #진동수 중심=0.05Hz

kx,ky=np.meshgrid(k,k)
fig=plt.figure()
ax=fig.subplots()
ax.pcolormesh(kx,ky,1/fks,cmap="Reds")
ax.scatter(x=p*np.cos(baz),y=p*np.sin(baz),s=200,c="blue",marker="+")
ax.set_xlabel(xlabel="kx(1/km)")
ax.set_ylabel(ylabel="ky(1/km)")
ax.grid(visible=True)
plt.show()
