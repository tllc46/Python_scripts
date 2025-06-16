#usage
#python xcorr.py 2020

import sys
from os.path import isfile

import numpy as np
import numpy.ma as ma
from scipy.signal import butter,ShortTimeFFT,detrend,sosfilt,hilbert,savgol_filter
from scipy.signal.windows import hann
import pandas as pd

from obspy import UTCDateTime,Stream,read

sec_day=86400 #[s]
epsilon=1e-10

sampling_rate=100 #[Hz]
delta=1/sampling_rate

#sub window
len_sub=48 #[s]
shift_sub=len_sub//2
npts_sub=len_sub*sampling_rate
npts_shift_sub=shift_sub*sampling_rate

#average window
nsub_avg=100
len_avg=shift_sub*nsub_avg
shift_avg=len_avg//2
npts_avg=(len_avg+len_sub//2)*sampling_rate

#day
year_ndays={2013:184,2014:365,2015:365,2016:366,2017:365,2018:365,2019:365,2020:366,2021:365,2022:365,2023:365}
offset_avg=3*3600 #[s]
year_start=int(sys.argv[1])
if year_start==2013:
    month_start=7
else:
    month_start=1
day_start=1
nday=year_ndays[year_start]
navg_day=1
day_udt=UTCDateTime(year=year_start,month=month_start,day=day_start)
navg=navg_day*nday

#frequency
nfreq=npts_sub//2+1
freqmin_filt=0.5 #[Hz]
freqmax_filt=1 #[Hz]

#stations
info_sta="center"
df=pd.read_csv(filepath_or_buffer=info_sta,sep=" ",names=["stnm","stla","stlo","stel"])
nsta=len(df)
ntriu=nsta*(nsta-1)//2
stnm_idx=dict(zip(df["stnm"],range(nsta)))
idx_triu=np.triu_indices(n=nsta,k=1)

#saving directory
path_save=f"/home/tllc46/48NAS1/tllc46/Aso/loc/{year_start}.npz"

#Butterworth bandpass filter
sos_butter=butter(N=4,Wn=[freqmin_filt,freqmax_filt],btype="band",fs=sampling_rate,output="sos")

#average window frequency normalization (whiten)
method_avg_fn="onebit"
win_dur=40 #[s]
len_fn_savgol=11
if method_avg_fn:
    win_size=win_dur*sampling_rate
    win=hann(M=win_size,sym=False)
    hop=win_size//2
    stf_avg=ShortTimeFFT(win=win,hop=hop,fs=sampling_rate,phase_shift=None)
    p0,p1=stf_avg.p_range(n=npts_avg)
    spctr_avg=np.empty(shape=(nsta,win_size//2+1,p1-p0),dtype=complex)

#average window time normalization
method_avg_tn="smooth"
len_tn_savgol=901

#sub window frequency normalization (whiten)
method_sub_fn=None
win=hann(M=npts_sub,sym=False)
stf_sub=ShortTimeFFT(win=win,hop=npts_shift_sub,fs=sampling_rate,phase_shift=None)
p0=stf_sub.lower_border_end[1]
p1=stf_sub.upper_border_begin(n=npts_avg)[1]

#arrays
stnm_avg=np.empty(shape=nsta,dtype=bool)
data_avg=np.empty(shape=(nsta,npts_avg))
spctr_sub=np.empty(shape=(nsta,nfreq,nsub_avg),dtype=complex)
cov_avg=np.empty(shape=(nsta,nsta,nfreq),dtype=complex)
cov_avg_diag=np.empty(shape=(nsta,nfreq))
denom=np.empty(shape=nsta)
xcorr=np.empty(shape=(navg,ntriu,npts_sub))

#stream
st=Stream()

def bandpass():
    data_avg[:nsta_avg,:]=sosfilt(sos=sos_butter,x=data_avg[:nsta_avg])
    data_avg[:nsta_avg,:]=sosfilt(sos=sos_butter,x=data_avg[:nsta_avg,::-1])[:nsta_avg,::-1]

def avg_fn():
    if not method_avg_fn:
        return

    spctr_avg[:nsta_avg,:,:]=stf_avg.stft(x=data_avg[:nsta_avg])
    if method_avg_fn=="onebit":
        spctr_avg[:nsta_avg,:,:]/=abs(spctr_avg[:nsta_avg])+epsilon
    elif method_avg_fn=="smooth":
        spctr_avg[:nsta_avg,:,:]/=savgol_filter(x=abs(spctr_avg[:nsta_avg]),window_length=len_fn_savgol,polyorder=1,axis=1)+epsilon
    data_avg[:nsta_avg,:]=stf_avg.istft(S=spctr_avg[:nsta_avg])[:,:npts_avg]

def avg_tn():
    if method_avg_tn=="onebit":
        data_avg[:nsta_avg,:]=np.sign(data_avg[:nsta_avg])
        ### original covseisnet ###
        #data_avg[:nsta_avg,:]/=abs(data_avg[:nsta_avg])+epsilon
    elif method_avg_tn=="smooth":
        data_avg[:nsta_avg,:]/=savgol_filter(x=abs(hilbert(x=data_avg[:nsta_avg])),window_length=len_tn_savgol,polyorder=1)+epsilon

def pre_proc():
    data_avg[:nsta_avg,:]=detrend(data=data_avg[:nsta_avg])
    bandpass()
    avg_fn()
    avg_tn()

def sub_fn():
    if method_sub_fn=="sub":
        spctr_sub[:nsta_avg,:,:]/=abs(spctr_sub[:nsta_avg])+epsilon
    elif method_sub_fn=="avg":
        spctr_sub[:nsta_avg,:,:]/=np.mean(a=abs(spctr_sub[:nsta_avg]),axis=2)[:,:,None]

def calc_cov():
    spctr_sub[:nsta_avg,:,:]=stf_sub.stft_detrend(x=data_avg[:nsta_avg],detr="linear",p0=p0,p1=p1)
    sub_fn()
    cov_avg[:nsta_avg,:nsta_avg,:]=np.einsum("ift,jft->ijf",spctr_sub[:nsta_avg],np.conj(spctr_sub[:nsta_avg]))/nsub_avg
    cov_avg_diag[:nsta_avg,:]=np.real(val=np.diagonal(a=cov_avg[:nsta_avg,:nsta_avg]).T)
    denom[:nsta_avg]=(2*np.sum(a=cov_avg_diag[:nsta_avg,1:-1],axis=1)+cov_avg_diag[:nsta_avg,0]+cov_avg_diag[:nsta_avg,-1])/npts_sub
    cov_avg[:nsta_avg,:nsta_avg,:]/=np.sqrt(denom[:nsta_avg,None,None]*denom[:nsta_avg,None])

def calc_xcorr():
    xcorr[idx_avg,idx_flat_avg,:]=np.fft.fftshift(x=np.fft.irfft(a=cov_avg[idx_triu_avg]),axes=1)
    xcorr[idx_avg,~idx_flat_avg,:]=np.nan

def calc_day():
    global nsta_avg,idx_avg,idx_flat_avg,idx_triu_avg

    for i in range(navg_day):
        nsta_avg=0
        stnm_avg[:]=0
        idx_avg=idx_day*navg_day+i
        start_avg=day_udt+offset_avg+i*shift_avg
        st_avg=st.slice(starttime=start_avg,endtime=start_avg+len_avg+len_sub//2-delta)
        nsta_avg_0=len(st_avg)

        for j in range(nsta_avg_0):
            tr=st_avg[j]
            if len(tr.data)<npts_avg or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_avg:
                continue

            stnm_avg[stnm_idx[tr.stats.station+"."+tr.stats.channel]]=True
            data_avg[nsta_avg,:]=tr.data
            nsta_avg+=1

        idx_flat_avg=(stnm_avg[:,None] & stnm_avg)[idx_triu]
        idx_triu_avg=np.triu_indices(n=nsta_avg,k=1)

        pre_proc()
        calc_cov()
        calc_xcorr()

def main():
    global idx_day,st,day_udt

    for idx_day in range(nday):
        for j in range(nsta):
            stnm=df.loc[j,"stnm"]
            path="/home/tllc46/48NAS1/tllc46/Aso_data/seismograms/"+day_udt.strftime(format="%Y.%j")+"/"+stnm+".SAC"
            if isfile(path=path):
                st+=read(pathname_or_url=path,format="SAC",byteorder="little")

        calc_day()
        st.clear()
        print(day_udt,"done")
        day_udt+=sec_day

    np.savez(file=path_save,xcorr=xcorr)

main()
