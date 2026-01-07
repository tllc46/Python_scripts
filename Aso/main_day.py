#usage
#python main_day.py t04 x195 2014-11

import sys
from os import makedirs
from os.path import isdir,isfile
from importlib import import_module
from datetime import datetime

import numpy as np
import numpy.ma as ma
from numpy.linalg import eigvalsh,eigh
from scipy.signal import cheb2ord,cheby2,sosfilt,butter,ShortTimeFFT,savgol_filter,hilbert,detrend
from scipy.signal.windows import hann
import pandas as pd

from obspy import Stream,read

sys.path.append("/home/tllc46/Aso/tstep")

mdl_t=import_module(name=sys.argv[1])
mdl_p=import_module(name=sys.argv[2])

calc_type="xcorr"
method_cohrnc=None
method_xcorr="default"
bool_deci=False
bool_fltr=True

#constant
sec_day=86400 #[s]
epsilon=1e-10

#sampling rate
sampling_rate_0=100 #[Hz]
delta=1/sampling_rate_0
if bool_deci:
    sampling_rate=25 #[Hz]
    factor=sampling_rate_0//sampling_rate
else:
    sampling_rate=sampling_rate_0 #[Hz]

#sub window
len_sub=mdl_p.len_sub
shift_sub=len_sub//2
npts_shift_sub=shift_sub*sampling_rate
ovrlp_sub=shift_sub
npts_sub=len_sub*sampling_rate

#average window
nsub_avg=mdl_p.nsub_avg
len_avg=shift_sub*nsub_avg+ovrlp_sub
npts_avg=len_avg*sampling_rate
npts_avg_0=len_avg*sampling_rate_0
shift_avg=mdl_t.shift_avg
offset_avg=mdl_t.offset_avg
navg_day=mdl_t.navg_day
if sec_day<offset_avg+shift_avg*(navg_day-1)+len_avg:
    st_step="conti"
    ovrlp_avg=(offset_avg+shift_avg*(navg_day-1)+len_avg)-sec_day
else:
    st_step="day"

#day
mdl_t.init(str_date=sys.argv[3])
udt_cur=mdl_t.udt_cur
nday=mdl_t.nday
udt_next=udt_cur+sec_day
navg=mdl_t.navg
npts_day=(sec_day+ovrlp_avg)*sampling_rate_0

#frequency
if calc_type=="cohrnc":
    idx_fmin=int(np.floor(mdl_p.fmin*len_sub))
    idx_fmax=int(np.ceil(mdl_p.fmax*len_sub))+1
    nfreq=idx_fmax-idx_fmin
elif calc_type=="xcorr":
    nfreq=npts_sub//2+1

#station
df=pd.read_csv(filepath_or_buffer=mdl_p.info_sta,sep=" ",names=["stnm","stla","stlo","stel"])
nsta=len(df)
ntriu=nsta*(nsta-1)//2
idx_triu=np.triu_indices(n=nsta,k=1)
stnm_idx=dict(zip(df["stnm"],range(nsta)))

#decimate
def init_deci():
    global sos_cheb

    order=1e99
    wp=1/factor
    while True:
        if order<=12:
            break
        wp*=0.99
        order,wn=cheb2ord(wp=wp,ws=1/factor,gpass=1,gstop=96)
    sos_cheb=cheby2(N=order,rs=96,Wn=wn,output="sos")
if bool_deci:
    init_deci()
def decimate():
    data_avg_0[:nsta_day,:]=sosfilt(sos=sos_cheb,x=data_avg_0[:nsta_day])
    data_avg[:nsta_day,:]=data_avg_0[:nsta_day,::factor]

#bandpass filter
if bool_fltr:
    sos_butter=butter(N=4,Wn=[mdl_p.fmin,mdl_p.fmax],btype="band",fs=sampling_rate,output="sos")
def bandpass():
    data_avg[:nsta_day,:]=sosfilt(sos=sos_butter,x=data_avg[:nsta_day])
    data_avg[:nsta_day,:]=sosfilt(sos=sos_butter,x=data_avg[:nsta_day,::-1])[:,::-1]

#average window frequency normalization (whiten)
method_avg_fn=mdl_p.method_avg_fn
len_fn_savgol=mdl_p.len_fn_savgol
if method_avg_fn:
    win_size=mdl_p.win_dur_fn*sampling_rate
    win=hann(M=win_size,sym=False)
    hop=win_size//2
    stf_avg=ShortTimeFFT(win=win,hop=hop,fs=sampling_rate,phase_shift=None)
    p0,p1=stf_avg.p_range(n=npts_avg)
    spctr_avg=np.empty(shape=(nsta,win_size//2+1,p1-p0),dtype=complex)
def avg_fn():
    if not method_avg_fn:
        return

    spctr_avg[:nsta_day,:,:]=stf_avg.stft(x=data_avg[:nsta_day])
    if method_avg_fn=="onebit":
        spctr_avg[:nsta_day,:,:]/=abs(spctr_avg[:nsta_day])+epsilon
    elif method_avg_fn=="smooth":
        spctr_avg[:nsta_day,:,:]/=savgol_filter(x=abs(spctr_avg[:nsta_day]),window_length=len_fn_savgol,polyorder=1,axis=1)+epsilon
    data_avg[:nsta_day,:]=stf_avg.istft(S=spctr_avg[:nsta_day])[:,:npts_avg]

#average window time normalization (whiten)
method_avg_tn=mdl_p.method_avg_tn
len_tn_savgol=mdl_p.len_tn_savgol
def avg_tn():
    if method_avg_tn=="onebit":
        data_avg[:nsta_day,:]=np.sign(data_avg[:nsta_day])
        ### original covseisnet ###
        #data_avg[:nsta_day,:]/=abs(data_avg[:nsta_day])+epsilon
    elif method_avg_tn=="smooth":
        data_avg[:nsta_day,:]/=savgol_filter(x=abs(hilbert(x=data_avg[:nsta_day])),window_length=len_tn_savgol,polyorder=1)+epsilon

#sub window frequency normalization (whiten)
method_sub_fn=mdl_p.method_sub_fn
win=hann(M=npts_sub,sym=False)
stf_sub=ShortTimeFFT(win=win,hop=npts_shift_sub,fs=sampling_rate,phase_shift=None)
p0=stf_sub.lower_border_end[1]
p1=stf_sub.upper_border_begin(n=npts_avg)[1]
def sub_fn():
    if calc_type=="cohrnc":
        spctr_sub[:nsta_day,:,:]=stf_sub.stft_detrend(x=data_avg[:nsta_day],detr="linear",p0=p0,p1=p1)[:,idx_fmin:idx_fmax]
    elif calc_type=="xcorr":
        spctr_sub[:nsta_day,:,:]=stf_sub.stft_detrend(x=data_avg[:nsta_day],detr="linear",p0=p0,p1=p1)

    if method_sub_fn=="sub":
        spctr_sub[:nsta_day,:,:]/=abs(spctr_sub[:nsta_day])+epsilon
    elif method_sub_fn=="avg":
        spctr_sub[:nsta_day,:,:]/=np.mean(a=abs(spctr_sub[:nsta_day]),axis=2,keepdims=True)

#stream
st_cur=Stream()
st_next=Stream()
def read_cur():
    global st_cur

    for i in range(nsta):
        stnm=df.loc[i,"stnm"]
        path="/home/tllc46/48NAS1/tllc46/Aso_data/"+udt_cur.strftime(format="%Y.%j")+"/"+stnm+".SAC"
        if isfile(path=path):
            st_cur+=read(pathname_or_url=path,format="SAC",byteorder="little")
def read_next():
    global st_next

    for i in range(nsta):
        stnm=df.loc[i,"stnm"]
        path="/home/tllc46/48NAS1/tllc46/Aso_data/"+udt_next.strftime(format="%Y.%j")+"/"+stnm+".SAC"
        if isfile(path=path):
            st_next+=read(pathname_or_url=path,format="SAC",byteorder="little")

#array
data_day=np.empty(shape=(nsta,npts_day))
data_avg=np.empty(shape=(nsta,npts_avg))
spctr_sub=np.empty(shape=(nsta,nfreq,nsub_avg),dtype=complex)

if calc_type=="cohrnc":
    cohrnc_day=np.empty(shape=(nfreq,navg_day))
    cohrnc=np.empty(shape=(nfreq,nday))
elif calc_type=="xcorr":
    stnm_day=np.empty(shape=nsta,dtype=bool)
    xcorr_day=np.empty(shape=(navg_day,ntriu,npts_sub))
    xcorr=np.empty(shape=(nday,ntriu,npts_sub))
    if st_step=="conti":
        order=np.empty(shape=nsta,dtype=int)
        idx_order=np.empty(shape=nsta,dtype=int)

if method_cohrnc=="simple":
    cov_sub=np.empty(shape=(nsta,nsta,nfreq,nsub_avg),dtype=complex)
elif method_cohrnc=="avg":
    cov_avg=np.empty(shape=(nfreq,nsta,nsta),dtype=complex)
    cov_avg_diag=np.empty(shape=(nfreq,nsta))
elif method_cohrnc=="eig":
    cov_avg=np.empty(shape=(nfreq,nsta,nsta),dtype=complex)
    eig_val=np.empty(shape=(nfreq,nsta))
elif method_xcorr=="default":
    cov_avg=np.empty(shape=(nfreq,nsta,nsta),dtype=complex)
    cov_avg_diag=np.empty(shape=(nfreq,nsta))
    denom=np.empty(shape=nsta)
elif method_xcorr=="eig":
    cov_avg=np.empty(shape=(nfreq,nsta,nsta),dtype=complex)
    eig_val=np.empty(shape=(nfreq,nsta))
    eig_vec=np.empty(shape=(nfreq,nsta,nsta),dtype=complex)

if bool_deci:
    data_avg_0=np.empty(shape=(nsta,npts_avg_0))

#log file
if not isdir("logs"):
    makedirs(name="logs")
file=open(file="logs/"+sys.argv[2]+"."+sys.argv[3],mode="w",buffering=1)
def print_status(idx_avg):
    now=datetime.now()
    now_str=now.strftime(format="[%m-%dT%H:%M:%S]")
    file.write(now_str+f" {idx_avg+1:03}/{navg} done\n")

#saving directory
if calc_type=="cohrnc":
    path_save="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/cohrnc"
elif calc_type=="xcorr":
    path_save="/home/tllc46/48NAS1/tllc46/Aso/"+sys.argv[1]+"/xcorr"
if not isdir(path_save):
    makedirs(name=path_save)
path_save+="/"+sys.argv[2]+"."+sys.argv[3]+".npz"
if isfile(path=path_save):
    print("data already exists")
    exit()

def pre_proc():
    if bool_deci:
        data_avg_0[:nsta_day,:]=detrend(data=data_avg_0[:nsta_day])
        decimate()
    else:
        data_avg[:nsta_day,:]=detrend(data=data_avg[:nsta_day])
    if bool_fltr:
        bandpass()
    avg_fn()
    avg_tn()

def calc_cov():
    sub_fn()
    if method_cohrnc=="simple":
        cov_sub[:nsta_day,:nsta_day,:,:]=spctr_sub[:nsta_day,None]*np.conj(spctr_sub[:nsta_day])
    else:
        cov_avg[:,:nsta_day,:nsta_day]=np.einsum("ift,jft->fij",spctr_sub[:nsta_day],np.conj(spctr_sub[:nsta_day]))/nsub_avg

def calc_cohrnc(idx_avg_day):
    if method_cohrnc=="simple":
        cov_sub[:nsta_day,:nsta_day,:,:]=np.exp(1j*np.angle(z=cov_sub[:nsta_day,:nsta_day]))
        cohrnc_day[:,idx_avg_day]=np.mean(a=abs(np.mean(a=cov_sub[idx_triu_day],axis=2)),axis=0)
    elif method_cohrnc=="avg":
        cov_avg_diag[:,:nsta_day]=np.real(val=np.diagonal(a=cov_avg[:,:nsta_day,:nsta_day],axis1=1,axis2=2))
        cov_avg[:,:nsta_day,:nsta_day]/=np.sqrt(cov_avg_diag[:,:nsta_day,None]*cov_avg_diag[:,None,:nsta_day])
        cohrnc_day[:,idx_avg_day]=np.nanmean(a=abs(cov_avg[:,idx_triu_day[0],idx_triu_day[1]]),axis=1) #"nan"mean: gap은 아니지만, 0만으로 이루어진 data가 있다면, 0으로 나누어져 nan 발생
    elif method_cohrnc=="eig":
        eig_val[:,:nsta_day]=eigvalsh(a=cov_avg[:,:nsta_day,:nsta_day])
        cohrnc_day[:,idx_avg_day]=np.sum(a=np.arange(start=nsta_day-1,stop=-1,step=-1)*eig_val[:,:nsta_day],axis=1)
        cohrnc_day[:,idx_avg_day]/=np.sum(a=eig_val[:,:nsta_day],axis=1)

def calc_xcorr(idx_avg_day):
    if method_xcorr=="default":
        cov_avg_diag[:,:nsta_day]=np.real(val=np.diagonal(a=cov_avg[:,:nsta_day,:nsta_day],axis1=1,axis2=2))
        denom[:nsta_day]=2*np.sum(a=cov_avg_diag[1:-1,:nsta_day],axis=0)
        denom[:nsta_day]+=cov_avg_diag[0,:nsta_day]+cov_avg_diag[-1,:nsta_day]
        denom[:nsta_day]/=npts_sub
        cov_avg[:,:nsta_day,:nsta_day]/=np.sqrt(denom[:nsta_day,None]*denom[:nsta_day])
    elif method_xcorr=="eig":
        eig_val[:,:nsta_day],eig_vec[:,:nsta_day,:nsta_day]=eigh(a=cov_avg[:,:nsta_day,:nsta_day])
        cov_avg[:,:nsta_day,:nsta_day]=eig_vec[:,:nsta_day,nsta_day-1,None]*np.conj(eig_vec[:,None,:nsta_day,nsta_day-1])
        cov_avg[:,:nsta_day,:nsta_day]*=eig_val[:,nsta_day-1,None,None]

    xcorr_day[idx_avg_day,idx_flat_day,:]=np.fft.fftshift(x=np.fft.irfft(a=cov_avg[:,idx_triu_day[0],idx_triu_day[1]],axis=0),axes=0).T

def st2data():
    global nsta_day

    nsta_day_0=len(st_cur)
    j=0
    for i in range(nsta_day_0):
        tr=st_cur[i]
        if len(tr.data)<npts_day or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_day:
            continue
        data_day[j,:]=tr.data

        if calc_type=="xcorr":
            stnm_day[stnm_idx[tr.stats.station+"."+tr.stats.channel]]=True

        j+=1

    nsta_day=j

def st2data_order():
    global nsta_day

    i=0
    while i<len(st_cur):
        tr=st_cur[i]
        if len(tr.data)<npts_day or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_day:
            st_cur.pop(index=i)
            continue

        stnm_day[stnm_idx[tr.stats.station+"."+tr.stats.channel]]=True
        order[i]=stnm_idx[tr.stats.station+"."+tr.stats.channel]
        i+=1

    nsta_day=len(st_cur)
    idx_order[:nsta_day]=np.argsort(a=np.argsort(a=order[:nsta_day]))
    for i in range(nsta_day):
        tr=st_cur[i]
        data_day[idx_order[i],:]=tr.data

def calc_day(idx_day):
    for idx_avg_day in range(navg_day):
        idx_avg=idx_day*navg_day+idx_avg_day
        idx_data=(offset_avg+idx_avg_day*shift_avg)*sampling_rate_0

        if bool_deci:
            data_avg_0[:nsta_day,:]=data_day[:nsta_day,idx_data:idx_data+npts_avg_0]
        else:
            data_avg[:nsta_day,:]=data_day[:nsta_day,idx_data:idx_data+npts_avg]

        pre_proc()
        calc_cov()
        if calc_type=="cohrnc":
            calc_cohrnc(idx_avg_day=idx_avg_day)
        elif calc_type=="xcorr":
            calc_xcorr(idx_avg_day=idx_avg_day)
        print_status(idx_avg=idx_avg)

    if calc_type=="cohrnc":
        cohrnc[:,idx_day]=np.nanmean(a=cohrnc_day,axis=1)
    elif calc_type=="xcorr":
        xcorr[idx_day,idx_flat_day,:]=np.nanmean(a=xcorr_day[:,idx_flat_day],axis=0)
        xcorr[idx_day,~idx_flat_day,:]=np.nan

def main_day():
    global udt_cur

    for idx_day in range(nday):
        read_cur()

        calc_day(idx_day=idx_day)
        st_cur.traces.clear()
        udt_cur+=sec_day

    if calc_type=="cohrnc":
        np.savez(file=path_save,cohrnc=cohrnc)
    elif calc_type=="xcorr":
        np.savez(file=path_save,xcorr=xcorr)
    file.close()

def main_conti():
    global st_cur,idx_triu_day,idx_flat_day,udt_cur,udt_next

    read_cur()
    for idx_day in range(nday):
        read_next()
        st_cur+=st_next.slice(endtime=udt_next+ovrlp_avg-delta)
        st_cur.merge(method=1)

        if calc_type=="xcorr":
            stnm_day[:]=0

        if st_step=="conti" and calc_type=="xcorr":
            st2data_order()
        else:
            st2data()

        idx_triu_day=np.triu_indices(n=nsta_day,k=1)
        if calc_type=="xcorr":
            idx_flat_day=(stnm_day[:,None] & stnm_day)[idx_triu]

        calc_day(idx_day=idx_day)
        st_cur.traces.clear()
        st_cur+=st_next
        st_next.traces.clear()
        udt_cur=udt_next
        udt_next+=sec_day

    if calc_type=="cohrnc":
        np.savez(file=path_save,cohrnc=cohrnc)
    elif calc_type=="xcorr":
        np.savez(file=path_save,xcorr=xcorr)
    file.close()

if st_step=="day":
    main_day()
elif st_step=="conti":
    main_conti()
