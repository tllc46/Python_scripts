#implementation of scipy.signal.ShortTimeFFT

import numpy as np
from scipy.signal.windows import hann

from obspy import read

st=read(pathname_or_url="05.SS05.HHZ.2013.304",format="MSEED",header_byteorder="<")
data=st[0].data[:360000]

win=hann(M=40)
hop=20

#ShortTimeFFT
#equivalent to ShortTimeFFT(win=win,hop=hop,phase_shift=None)
n=len(data)

m_num=len(win)
m_num_mid=m_num//2
nonzl=np.flatnonzero(a=win)[0]
nonzr=(m_num-1)-np.flatnonzero(a=win)[-1]
win_l=-m_num_mid
win_zl=win_l+nonzl #first non-zero window element index
win_r=(m_num-1)-m_num_mid
win_zr=win_r-nonzr #last non-zero window element index

w2=win**2
DD=w2.copy()
for k in range(hop,m_num,hop):
    DD[k:]+=w2[:-k]
    DD[:-k]+=w2[k:]
dual_win=win/DD

pmin=-(win_r//hop) #0<=win_r+p*hop
plb=np.ceil((-win_zl)/hop) #0<=win_zl+p*hop
pmax=np.ceil((n-win_zl)/hop) #n<=win_zl+p*hop
pub=np.ceil((n-win_zr)/hop) #n<=win_zr+p*hop
kmin=pmin*hop-m_num_mid

#ShortTimeFFT.stft()
#equivalent to S=ShortTimeFFT.stft(x=data,p0=p0,p1=p1)
p0=pmin
p1=pmax

k0=p0*hop-m_num_mid
k1=(p1-1)*hop+(m_num-m_num_mid)
x1=np.zeros(shape=k1-k0)

left_diff=0-k0
right_diff=n-k1
npts=n

if 0<=left_diff:
    i0=0
    j0=left_diff
else:
    i0=-left_diff
    j0=0
    npts+=left_diff

if 0<=right_diff:
    npts-=right_diff

x1[j0:j0+npts]=data[i0:i0+npts]

f_pts=m_num//2+1
S=np.empty(shape=(f_pts,p1-p0),dtype=complex)
for i in range(p1-p0):
    segment=x1[i*hop:i*hop+m_num]*win
    S[:,i]=np.fft.rfft(a=segment)

#ShortTimeFFT.istft()
#equivalent to x=ShortTimeFFT.istft(S=S)
x=np.zeros(shape=k1-k0)
for i in range(p1-p0):
    xs=np.fft.irfft(a=S[:,i],n=m_num)*dual_win
    x[i*hop:i*hop+m_num]+=xs
x=x[-kmin:]
