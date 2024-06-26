from obspy import read
from obspy.core import AttribDict
from obspy.signal.array_analysis import get_geometry
from obspy.signal.util import next_pow_2
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

import numpy as np
import matplotlib.pyplot as plt

d2km=180/(np.pi*6371)
sll_x=-3 #최소 x slowness
slm_x=5 #최대 x slowness
sll_y=10 #최소 y slowness
slm_y=20 #최대 y slowness
sl_s=0.1 #slowness 단계
frqlow=1 #최소 진동수
frqhigh=2 #최대 진동수
frqfoc=1.5 #focusing 진동수
ssm="incoherent" #신호 부공간 방법, "coherent" or "incoherent"
method="Bartlett" #beamforming 방법, "Bartlett" or "Capon" or "MUSIC"

def fk_cssm(sl_x,sl_y):
    steer=np.exp(1j*2*np.pi*frqfoc*d2km*(sl_x*geometry[:,0]+sl_y*geometry[:,1]))
    tmp=np.matmul(scm,steer)
    power=np.real(val=np.matmul(steer.conj(),tmp))
    return power

def fk_issm(sl_x,sl_y):
    power=0
    for i in range(nf):
        steer=np.exp(1j*2*np.pi*frq[nlow+i]*d2km*(sl_x*geometry[:,0]+sl_y*geometry[:,1]))
        tmp=np.matmul(scm[:,:,i],steer)
        tmp2=abs(np.matmul(steer.conj(),tmp))
        if method=="Capon":
            tmp2=1/tmp2
        power+=tmp2
    return power

np_fk_cssm=np.frompyfunc(fk_cssm,nin=2,nout=1)
np_fk_issm=np.frompyfunc(fk_issm,nin=2,nout=1)

#main
#1. data 읽기
st=read(pathname_or_url="cutdata/24/Z/*",format="SAC",byteorder="little")
nch=len(st)

#2. 지진계 중심으로부터 거리 vector
for i in range(nch):
    st[i].stats.coordinates=AttribDict({"latitude":st[i].stats.sac.stla,"longitude":st[i].stats.sac.stlo,"elevation":-12345})
geometry=get_geometry(stream=st,return_center=True)
center=geometry[-1]
geometry=geometry[:-1]
_,_,baz=gps2dist_azimuth(lat1=st[0].stats.sac.evla,lon1=st[0].stats.sac.evlo,lat2=center[1],lon2=center[0])
angle=0.5*np.pi-np.radians(baz)

#3. 참조 slowness 계산
model=TauPyModel(model="ak135")
arrivals=model.get_travel_times_geo(source_depth_in_km=st[0].stats.sac.evdp,source_latitude_in_deg=st[0].stats.sac.evla,source_longitude_in_deg=st[0].stats.sac.evlo,receiver_latitude_in_deg=center[1],receiver_longitude_in_deg=center[0],phase_list=["P","Pn"])
sl_P=arrivals[2].ray_param_sec_degree
sl_Pn=arrivals[1].ray_param_sec_degree
phase=sl_P

#4. rfft 준비
nsamples=st[0].stats.npts
fs=st[0].stats.sampling_rate
nfft=next_pow_2(i=nsamples)
deltaf=fs/nfft
nlow=int(frqlow/deltaf+0.5)
nhigh=int(frqhigh/deltaf+0.5)
nlow=max(1,nlow) #최소 진동수가 0Hz보다는 커야 한다
nhigh=min(nfft//2-1,nhigh) #최대 진동수가 Nyquist 진동수보다는 작아야 한다
nf=nhigh-nlow+1 #최소 진동수와 최대 진동수도 포함
frq=np.fft.rfftfreq(n=nfft//2+1,d=1/fs) #진동수 범위

#5. rfft
output=np.empty(shape=(nch,nf),dtype=np.complex128)
for i in range(nch):
    dat=st[i].data
    output[i]=np.fft.rfft(a=dat,n=nfft)[nlow:nlow+nf]

#6. 공분산 행렬 계산
if ssm=="coherent":
    #6-1. focus
    for i in range(nf):
        output[:,i]*=np.exp(1j*2*np.pi*(frqfoc-frq[nlow+i])*d2km*phase*(np.cos(angle)*geometry[:,0]+np.sin(angle)*geometry[:,1]))

    #6-2. 공분산 행렬의 진동수에 대한 평균 계산
    scm=np.empty(shape=(nch,nch),dtype=np.complex128)
    for i in range(nch):
        for j in range(i,nch):
            scm[i,j]=sum(output[i]*output[j].conj())
            if i!=j:
                scm[j,i]=scm[i,j].conj()
else:
    scm=np.empty(shape=(nch,nch,nf),dtype=np.complex128)
    for i in range(nch):
        for j in range(i,nch):
            scm[i,j]=output[i]*output[j].conj()
            if method=="Capon":
                scm[i,j]/=np.abs(scm[i,j].sum())
            if i!=j:
                scm[j,i]=scm[i,j].conj()

#7. 공분산 행렬에 beamforming 방법 적용
if ssm=="coherent":
    if method=="Capon":
        w,z=np.linalg.eig(a=scm)
        tmp=np.matmul(z,np.diag(v=1/w))
        scm=np.matmul(tmp,np.conjugate(np.transpose(a=z)))
    
    elif method=="MUSIC":
        w,z=np.linalg.eig(a=scm)
        ind=np.argsort(a=np.real(val=w))
        ind=ind[:-1]
        scm=np.matmul(z[:,ind],np.conjugate(np.transpose(a=z[:,ind])))
else:
    if method=="Capon":
    #P(f)=1/(e.H R(f)^-1 e)
        for i in range(nf):
           scm[:,:,i]=np.linalg.pinv(a=scm[:,:,i],rcond=1e-6)

#8. beamforming
sl_x=np.arange(stop=(slm_x-sll_x)/sl_s)*sl_s+sll_x
sl_y=np.arange(stop=(slm_y-sll_y)/sl_s)*sl_s+sll_y
sl_x,sl_y=np.meshgrid(sl_x,sl_y)
if ssm=="coherent":
    fks=np_fk_cssm(sl_x=sl_x,sl_y=sl_y)
    fks=np.array(object=fks,dtype=float)
    if method=="Capon" or method=="MUSIC":
        fks=1/fks
else:
    fks=np_fk_issm(sl_x=sl_x,sl_y=sl_y)
    fks=np.array(object=fks,dtype=float)

#9. find maximum index
ix,iy=np.unravel_index(indices=fks.argmax(),shape=fks.shape)

#10. plot
plt.figure()
ax=fig.subplots()
ax.pcolormesh(sl_x,sl_y,fks,cmap="Reds")
ax.scatter(x=sl_x[ix,iy],y=sl_y[ix,iy],s=200,c="green",marker="+")
ax.scatter(x=sl_Pn*np.cos(angle),y=sl_Pn*np.sin(angle),s=200,c="blue",marker="+")
ax.scatter(x=sl_P*np.cos(angle),y=sl_P*np.sin(angle),s=200,c="blue",marker="+")
ax.set_xlabel(xlabel="sl_x(s/°)")
ax.set_ylabel(ylabel="sl_y(s/°)")
ax.grid(visible=True)
plt.show()
