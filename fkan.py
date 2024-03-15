from obspy import read
from obspy.core import AttribDict
from obspy.signal.array_analysis import get_geometry
from obspy.signal.util import next_pow_2
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

import numpy as np
import matplotlib.pyplot as plt

d2km=180/(np.pi*6371)
sll_x=-3 #slowness least x
slm_x=5 #slowness max x
sll_y=10 #slowness least y
slm_y=20 #slowness max y
sl_s=0.1 #slowness step
frqlow=1 #frequency low
frqhigh=2 #frequency high
frqfoc=1.5 #frequency focusing
ssm="incoherent" #signal subspace method, "coherent" or "incoherent"
method="Bartlett" #beamforming method, "Bartlett" or "Capon" or "MUSIC"

def fk_cssm(sl_x,sl_y):
    steer=np.exp(1j*2*np.pi*frqfoc*d2km*(sl_x*geometry[:,0]+sl_y*geometry[:,1]))
    tmp=np.matmul(scm,steer)
    power=np.real(np.matmul(steer.conj(),tmp))
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

np_fk_cssm=np.frompyfunc(fk_cssm,2,1)
np_fk_issm=np.frompyfunc(fk_issm,2,1)

#main
#1. read data
st=read(pathname_or_url="cutdata/24/Z/*",format="SAC",byteorder="little")
#st=st[30:90]
nch=len(st)

#2. get array geometry
for i in range(nch):
    st[i].stats.coordinates=AttribDict({"latitude":st[i].stats.sac.stla,"longitude":st[i].stats.sac.stlo,"elevation":-12345})
geometry=get_geometry(stream=st,return_center=True)
center=geometry[-1]
geometry=geometry[:-1]
_,_,baz=gps2dist_azimuth(st[0].stats.sac.evla,st[0].stats.sac.evlo,center[1],center[0])
angle=0.5*np.pi-baz*np.pi/180

#3. get reference slowness
model=TauPyModel(model="ak135")
arrivals=model.get_travel_times_geo(st[0].stats.sac.evdp,st[0].stats.sac.evla,st[0].stats.sac.evlo,center[1],center[0],phase_list=["P","Pn"])
sl_P=arrivals[2].ray_param_sec_degree
sl_Pn=arrivals[1].ray_param_sec_degree
phase=sl_P

#5. generate plan for rfft
nsamples=st[0].stats.npts
fs=st[0].stats.sampling_rate
nfft=next_pow_2(nsamples)
deltaf=fs/nfft
nlow=int(frqlow/deltaf+0.5)
nhigh=int(frqhigh/deltaf+0.5)
nlow=max(1,nlow) #avoid using the offset
nhigh=min(nfft//2-1,nhigh) #avoid using Nyquist
nf=nhigh-nlow+1 #include upper and lower frequency
frq=np.linspace(0,fs/2,nfft//2+1) #frequency range

#6. rfft
output=np.empty((nch,nf),dtype=np.complex128)
for i in range(nch):
    dat=st[i].data
    output[i]=np.fft.rfft(dat,nfft)[nlow:nlow+nf]

#7. computing the covariances of the signal at different receivers
if ssm=="coherent":
    #7-1. focus
    for i in range(nf):
        output[:,i]*=np.exp(1j*2*np.pi*(frqfoc-frq[nlow+i])*d2km*phase*(np.cos(angle)*geometry[:,0]+np.sin(angle)*geometry[:,1]))

    #7-2. compute average covariance matrix
    scm=np.empty((nch,nch),dtype=np.complex128)
    for i in range(nch):
        for j in range(i,nch):
            scm[i,j]=sum(output[i]*output[j].conj())
            if i!=j:
                scm[j,i]=scm[i,j].conj()
else:
    scm=np.empty((nch,nch,nf),dtype=np.complex128)
    for i in range(nch):
        for j in range(i,nch):
            scm[i,j]=output[i]*output[j].conj()
            if method=="Capon":
                scm[i,j]/=np.abs(scm[i,j].sum())
            if i!=j:
                scm[j,i]=scm[i,j].conj()

#8. apply methods
if ssm=="coherent":
    if method=="Capon":
        w,z=np.linalg.eig(scm)
        tmp=np.matmul(z,np.diag(1/w))
        scm=np.matmul(tmp,np.conjugate(np.transpose(z)))
    
    elif method=="MUSIC":
        w,z=np.linalg.eig(scm)
        ind=np.argsort(np.real(w))
        ind=ind[:-1]
        scm=np.matmul(z[:,ind],np.conjugate(np.transpose(z[:,ind])))
else:
    if method=="Capon":
    #P(f)=1/(e.H R(f)^-1 e)
        for i in range(nf):
           scm[:,:,i]=np.linalg.pinv(scm[:,:,i],rcond=1e-6)

#9. beamforming
sl_x=np.arange(sll_x,slm_x,sl_s)
sl_y=np.arange(sll_y,slm_y,sl_s)
sl_x,sl_y=np.meshgrid(sl_x,sl_y)
if ssm=="coherent":
    fks=np_fk_cssm(sl_x,sl_y)
    fks=np.array(fks,dtype=float)
    if method=="Capon" or method=="MUSIC":
        fks=1/fks
else:
    fks=np_fk_issm(sl_x,sl_y)
    fks=np.array(fks,dtype=float)

#10. find maximum index
ix,iy=np.unravel_index(fks.argmax(),fks.shape)

#11. plot
plt.figure()
plt.set_cmap("Reds")
plt.pcolormesh(sl_x,sl_y,fks)
plt.scatter(sl_x[ix,iy],sl_y[ix,iy],c="green",s=200,marker="+")
plt.scatter(sl_Pn*np.cos(angle),sl_Pn*np.sin(angle),c="blue",s=200,marker="+")
plt.scatter(sl_P*np.cos(angle),sl_P*np.sin(angle),c="blue",s=200,marker="+")
plt.xlabel("sl_x(s/°)")
plt.ylabel("sl_y(s/°)")
plt.grid()
plt.show()
