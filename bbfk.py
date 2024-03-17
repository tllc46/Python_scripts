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
    trace=np.real(np.trace(scm))
    for i in range(nch):
        v=np.real(scm[i][i])
        scale=np.sqrt(trace/(v*nch))
        for j in range(nch):
            scm[i][j]*=scale
            scm[j][i]*=scale

def regularize():
    trace=np.real(np.trace(scm))
    offset=trace*eps/nch
    for i in range(nch):
        scm[i][i]+=offset

def mlm():
    #fks/xbbfk.c/eigenanal()
    w,z=np.linalg.eig(scm)
    tmp=np.matmul(z,np.diag(1/w))
    sinv=np.matmul(tmp,np.conjugate(np.transpose(z)))
    return sinv

def music():
    #fks/xbbfk.c/eigenanal()
    w,z=np.linalg.eig(scm)
    ind=np.argsort(np.real(w))
    ind=ind[:-1]
    sinv=np.matmul(z[:,ind],np.conjugate(np.transpose(z[:,ind])))
    return sinv

def fkevalr():
    #fks/xbbfk.c/fkevalr()
    wv=np.linspace(-2*np.pi*wavenumber,2*np.pi*wavenumber,ssq)
    for i in range(ssq):
        for j in range(ssq):
            tmp=np.matmul(np.exp(-1j*(wv[i]*geometry[:,0]+wv[j]*geometry[:,1])),scm)
            fks[i][j]=np.real(np.matmul(tmp,np.exp(1j*(wv[i]*geometry[:,0]+wv[j]*geometry[:,1]))))

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
st=st.slice(start_com,end_com)
nsamples=st[0].stats.npts

for i in range(nch):
    st[i].stats.coordinates=AttribDict({"latitude":st[i].stats.sac.stla,"longitude":st[i].stats.sac.stlo,"elevation":-12345})

geometry=get_geometry(stream=st,return_center=True)
center=geometry[-1]
geometry=geometry[:-1]
_,_,baz=gps2dist_azimuth(st[0].stats.sac.evla,st[0].stats.sac.evlo,center[1],center[0])
baz*=np.pi/180

output=np.empty((nch,100),dtype=complex)
scm=np.zeros((nch,nch),dtype=complex)
fks=np.empty((ssq,ssq))
qstates=np.zeros((nch,11))

covmat()
scm=music()
fkevalr()

model=TauPyModel(model="ak135")
arrivals=model.get_travel_times_geo(st[0].stats.sac.evdp,st[0].stats.sac.evla,st[0].stats.sac.evlo,center[1],center[0],phase_list=["P"])
p=arrivals[0].ray_param_sec_degree
p=p*0.05/111.19 #frequency peak=0.05Hz

k=np.linspace(-wavenumber,wavenumber,ssq)
kx,ky=np.meshgrid(k,k)
plt.set_cmap("Reds")
plt.pcolormesh(kx,ky,1/fks)
plt.scatter(p*np.cos(baz),p*np.sin(baz),c="blue",s=200,marker="+")
plt.xlabel("kx(1/km)")
plt.ylabel("ky(1/km)")
plt.grid()
plt.show()
