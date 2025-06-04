import numpy as np
from scipy.signal import detrend,butter,sosfilt,cheb2ord,cheby2
from scipy.signal.windows import hann

sampling_rate=100
npts=48*sampling_rate

def detrend_obspy(data):
    data[:]=detrend(data=data)

def taper(data,max_percentage):
    wlen=int(max_percentage*npts)
    if 2*wlen==npts:
        taper_sides=hann(M=2*wlen)[:wlen]
    else:
        taper_sides=hann(M=2*wlen+1)[:wlen]

    data[:wlen]*=taper_sides
    data[-wlen:]*=taper_sides[::-1]

def bandpass(data,freqmin,freqmax):
    sos=butter(N=4,Wn=[freqmin,freqmax],btype="band",fs=sampling_rate,output="sos")
    data[:]=sosfilt(sos=sos,x=data)
    data[:]=sosfilt(sos=sos,x=data[::-1])[::-1]

def lowpass_cheby_2(factor):
    order=1e99
    wp=1/factor
    while True:
        if order<=12:
            break
        wp*=0.99
        order,wn=cheb2ord(wp=wp,ws=1/factor,gpass=1,gstop=96)

    sos=cheby2(N=order,rs=96,Wn=wn,output="sos")
    data[:]=sosfilt(sos=sos,x=data)

def decimate(data,factor):
    lowpass_cheby_2(factor=factor)
    return data[::factor]

def resample(data,factor):
    lowpass_cheby_2(factor=factor)
    large_w=np.fft.ifftshift(x=hann(M=npts,sym=False))[:npts//2+1]
    x=np.fft.rfft(a=data)*large_w
    f=np.fft.rfftfreq(n=npts,d=1/sampling_rate)
    large_f=np.fft.rfftfreq(n=npts//factor,d=factor/sampling_rate)
    large_y=np.interp(x=large_f,xp=f,fp=x)
    return np.fft.irfft(a=large_y)/factor
