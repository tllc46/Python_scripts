from obspy import read
from obspy.signal.util import next_pow_2

import numpy as np
from numpy.polynomial.polynomial import polyfit,polyval

st=read(pathname_or_url=f"/media/tllc46/data01/Donghae/cutdata_dec_NS/2/*",format="SAC",byteorder="little")
nsta=len(st)
npts=st[0].stats.npts
delta=st[0].stats.delta
lag=delta*np.arange(start=-npts+1,stop=npts)

npts_2pow=next_pow_2(i=2*npts-1)
nfft=npts_2pow//2+1
fourier=np.empty(shape=(nsta,nfft),dtype=complex)
factor=np.empty(shape=nsta)
L=np.zeros(shape=(nsta,nsta))
C=np.identity(n=nsta)

for i in range(nsta):
    tr=st[i]
    fourier[i]=np.fft.rfft(a=tr.data,n=npts_2pow)
    factor[i]=np.sqrt(sum(tr.data**2))

fourier_conj=np.conjugate(fourier)

for i in range(nsta-1):
    for j in range(i+1,nsta):
        cc=fourier[i]*fourier_conj[j]
        corr=np.real(val=np.fft.irfft(a=cc))
        corr=np.concatenate((corr[-npts+1:],corr[:npts]))
        test_max=max(corr[1:-1])
        index=np.argmax(a=corr[1:-1])+1
        coeff=polyfit(x=lag[index-1:index+2],y=corr[index-1:index+2],deg=2)
        time=-coeff[1]/(2*coeff[2])
        corr_max=polyval(x=time,c=coeff)/(factor[i]*factor[j])
        L[i,j]=time #양수일 경우, i행에 대해 j열 신호가 더 일찍 도착
        L[j,i]=-time
        C[i,j]=corr_max
        C[j,i]=corr_max
