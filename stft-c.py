import numpy as np

from obspy import read

st=read("05.SS05.HHZ.2013.304",format="MSEED",header_byteorder="<")
data=st[0].data[:360000]

nperseg=40

#stft
fac=np.arange(nperseg+1)
win=0.5*(1-np.cos(2*np.pi*fac/(nperseg)))[:-1]
nfft=nperseg
noverlap=nperseg//2
nstep=nperseg-noverlap
nadd=-len(data)%nstep
pad_length=nperseg+len(data)+nadd
x=np.zeros(shape=pad_length)
x[nperseg//2:-nperseg//2-nadd]=data
scale=1/(0.5*nperseg)
nseg=(pad_length-nperseg)//nstep+1
result=np.empty(shape=(nfft//2+1,nseg),dtype=complex)
for i in range(nseg):
    segment=x[i*nstep:i*nstep+nperseg]*win
    result[:,i]=np.fft.rfft(segment)
result*=scale

#istft
x=np.zeros(shape=pad_length)
norm=np.zeros(shape=pad_length)
for i in range(nseg):
    xsubs=np.fft.irfft(result[:,i])
    xsubs*=0.5*nperseg
    x[i*nstep:i*nstep+nperseg]+=xsubs*win
    norm[i*nstep:i*nstep+nperseg]+=win**2

x=x[nperseg//2:-nperseg//2-nadd]
norm=norm[nperseg//2:-nperseg//2-nadd]
x/=norm
