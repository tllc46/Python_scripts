#implementation of scipy.signal.stft() and scipy.signal.istft()

import numpy as np

from obspy import read

st=read(pathname_or_url="05.SS05.HHZ.2013.304",format="MSEED",header_byteorder="<")
data=st[0].data[:360000]

nperseg_f=40

#stft
#same as: result=scipy.signal.stft(x=data,nperseg=nperseg_f)[2]
fac_f=np.arange(stop=nperseg_f+1)
win_f=0.5*(1-np.cos(2*np.pi*fac_f/nperseg_f))[:-1]
nfft=nperseg_f
noverlap_f=nperseg_f//2
nstep_f=nperseg_f-noverlap_f
nadd=-(len(data)-(nperseg_f%2))%nstep_f
pad_length=(nperseg_f//2)*2+len(data)+nadd
x=np.zeros(shape=pad_length)
x[nperseg_f//2:-(nperseg_f//2)-nadd]=data
scale=1/(0.5*nperseg_f) #integral of hann window
nseg=(pad_length-nperseg_f+1)//nstep_f+bool((pad_length-nperseg_f+1)%nstep_f)
result=np.empty(shape=(nfft//2+1,nseg),dtype=complex)
for i in range(nseg):
    segment=x[i*nstep_f:i*nstep_f+nperseg_f]*win_f
    result[:,i]=np.fft.rfft(a=segment)
result*=scale

#istft
#same as: x=scipy.signal.istft(Zxx=result)[1]
nperseg_b=(nfft//2)*2
noverlap_b=nperseg_b//2
step_b=nperseg_b-noverlap_b
fac_b=np.arange(stop=nperseg_b+1)
win_b=0.5*(1-np.cos(2*np.pi*fac_b/nperseg_b))[:-1]
outputlength=nperseg_b+(nseg-1)*step_b
x=np.zeros(shape=outputlength)
norm=np.zeros(shape=outputlength)
for i in range(nseg):
    xsubs=np.fft.irfft(a=result[:,i])
    xsubs*=0.5*nperseg_b #integral of hann window
    x[i*step_b:i*step_b+nperseg_b]+=xsubs*win_b
    norm[i*step_b:i*step_b+nperseg_b]+=win_b**2

x=x[nperseg_b//2:-(nperseg_b//2)]
norm=norm[nperseg_b//2:-(nperseg_b//2)]
x/=norm
