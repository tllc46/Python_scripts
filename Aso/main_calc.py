#usage
#python main_calc.py t01 x01 2014

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
    data_avg_0[:nsta_avg,:]=sosfilt(sos=sos_cheb,x=data_avg_0[:nsta_avg])
    data_avg[:nsta_avg,:]=data_avg_0[:nsta_avg,::factor]

#bandpass filter
if bool_fltr:
    sos_butter=butter(N=4,Wn=[mdl_p.fmin,mdl_p.fmax],btype="band",fs=sampling_rate,output="sos")
def bandpass():
    data_avg[:nsta_avg,:]=sosfilt(sos=sos_butter,x=data_avg[:nsta_avg])
    data_avg[:nsta_avg,:]=sosfilt(sos=sos_butter,x=data_avg[:nsta_avg,::-1])[:,::-1]

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

    spctr_avg[:nsta_avg,:,:]=stf_avg.stft(x=data_avg[:nsta_avg])
    if method_avg_fn=="onebit":
        spctr_avg[:nsta_avg,:,:]/=abs(spctr_avg[:nsta_avg])+epsilon
    elif method_avg_fn=="smooth":
        spctr_avg[:nsta_avg,:,:]/=savgol_filter(x=abs(spctr_avg[:nsta_avg]),window_length=len_fn_savgol,polyorder=1,axis=1)+epsilon
    data_avg[:nsta_avg,:]=stf_avg.istft(S=spctr_avg[:nsta_avg])[:,:npts_avg]

#average window time normalization (whiten)
method_avg_tn=mdl_p.method_avg_tn
len_tn_savgol=mdl_p.len_tn_savgol
def avg_tn():
    if method_avg_tn=="onebit":
        data_avg[:nsta_avg,:]=np.sign(data_avg[:nsta_avg])
        ### original covseisnet ###
        #data_avg[:nsta_avg,:]/=abs(data_avg[:nsta_avg])+epsilon
    elif method_avg_tn=="smooth":
        data_avg[:nsta_avg,:]/=savgol_filter(x=abs(hilbert(x=data_avg[:nsta_avg])),window_length=len_tn_savgol,polyorder=1)+epsilon

#sub window frequency normalization (whiten)
method_sub_fn=mdl_p.method_sub_fn
win=hann(M=npts_sub,sym=False)
stf_sub=ShortTimeFFT(win=win,hop=npts_shift_sub,fs=sampling_rate,phase_shift=None)
p0=stf_sub.lower_border_end[1]
p1=stf_sub.upper_border_begin(n=npts_avg)[1]
def sub_fn():
    if calc_type=="cohrnc":
        spctr_sub[:nsta_avg,:,:]=stf_sub.stft_detrend(x=data_avg[:nsta_avg],detr="linear",p0=p0,p1=p1)[:,idx_fmin:idx_fmax]
    elif calc_type=="xcorr":
        spctr_sub[:nsta_avg,:,:]=stf_sub.stft_detrend(x=data_avg[:nsta_avg],detr="linear",p0=p0,p1=p1)

    if method_sub_fn=="sub":
        spctr_sub[:nsta_avg,:,:]/=abs(spctr_sub[:nsta_avg])+epsilon
    elif method_sub_fn=="avg":
        spctr_sub[:nsta_avg,:,:]/=np.mean(a=abs(spctr_sub[:nsta_avg]),axis=2,keepdims=True)

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
data_avg=np.empty(shape=(nsta,npts_avg))
spctr_sub=np.empty(shape=(nsta,nfreq,nsub_avg),dtype=complex)

if calc_type=="cohrnc":
    cohrnc=np.empty(shape=(nfreq,navg))
elif calc_type=="xcorr":
    stnm_avg=np.empty(shape=nsta,dtype=bool)
    xcorr=np.empty(shape=(navg,ntriu,npts_sub))
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
        data_avg_0[:nsta_avg,:]=detrend(data=data_avg_0[:nsta_avg])
        decimate()
    else:
        data_avg[:nsta_avg,:]=detrend(data=data_avg[:nsta_avg])
    if bool_fltr:
        bandpass()
    avg_fn()
    avg_tn()

def calc_cov():
    sub_fn()
    if method_cohrnc=="simple":
        cov_sub[:nsta_avg,:nsta_avg,:,:]=spctr_sub[:nsta_avg,None]*np.conj(spctr_sub[:nsta_avg])
    else:
        cov_avg[:,:nsta_avg,:nsta_avg]=np.einsum("ift,jft->fij",spctr_sub[:nsta_avg],np.conj(spctr_sub[:nsta_avg]))/nsub_avg

def calc_cohrnc(idx_avg):
    if method_cohrnc=="simple":
        cov_sub[:nsta_avg,:nsta_avg,:,:]=np.exp(1j*np.angle(z=cov_sub[:nsta_avg,:nsta_avg]))
        cohrnc[:,idx_avg]=np.mean(a=abs(np.mean(a=cov_sub[idx_triu_avg],axis=2)),axis=0)
    elif method_cohrnc=="avg":
        cov_avg_diag[:,:nsta_avg]=np.real(val=np.diagonal(a=cov_avg[:,:nsta_avg,:nsta_avg],axis1=1,axis2=2))
        cov_avg[:,:nsta_avg,:nsta_avg]/=np.sqrt(cov_avg_diag[:,:nsta_avg,None]*cov_avg_diag[:,None,:nsta_avg])
        cohrnc[:,idx_avg]=np.nanmean(a=abs(cov_avg[:,idx_triu_avg[0],idx_triu_avg[1]]),axis=1)
    elif method_cohrnc=="eig":
        eig_val[:,:nsta_avg]=eigvalsh(a=cov_avg[:,:nsta_avg,:nsta_avg])
        cohrnc[:,idx_avg]=np.sum(a=np.arange(start=nsta_avg-1,stop=-1,step=-1)*eig_val[:,:nsta_avg],axis=1)
        cohrnc[:,idx_avg]/=np.sum(a=eig_val[:,:nsta_avg],axis=1)

def calc_xcorr(idx_avg):
    if method_xcorr=="default":
        cov_avg_diag[:,:nsta_avg]=np.real(val=np.diagonal(a=cov_avg[:,:nsta_avg,:nsta_avg],axis1=1,axis2=2))
        denom[:nsta_avg]=2*np.sum(a=cov_avg_diag[1:-1,:nsta_avg],axis=0)
        denom[:nsta_avg]+=cov_avg_diag[0,:nsta_avg]+cov_avg_diag[-1,:nsta_avg]
        denom[:nsta_avg]/=npts_sub
        cov_avg[:,:nsta_avg,:nsta_avg]/=np.sqrt(denom[:nsta_avg,None]*denom[:nsta_avg])
    elif method_xcorr=="eig":
        eig_val[:,:nsta_avg],eig_vec[:,:nsta_avg,:nsta_avg]=eigh(a=cov_avg[:,:nsta_avg,:nsta_avg])
        cov_avg[:,:nsta_avg,:nsta_avg]=eig_vec[:,:nsta_avg,nsta_avg-1,None]*np.conj(eig_vec[:,None,:nsta_avg,nsta_avg-1])
        cov_avg[:,:nsta_avg,:nsta_avg]*=eig_val[:,nsta_avg-1,None,None]

    xcorr[idx_avg,idx_flat_avg,:]=np.fft.fftshift(x=np.fft.irfft(a=cov_avg[:,idx_triu_avg[0],idx_triu_avg[1]],axis=0),axes=0).T
    xcorr[idx_avg,~idx_flat_avg,:]=np.nan

def st2data():
    global nsta_avg

    nsta_avg_0=len(st_avg)
    j=0
    for i in range(nsta_avg_0):
        tr=st_avg[i]
        if bool_deci:
            if len(tr.data)<npts_avg_0 or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_avg_0:
                continue
            data_avg_0[j,:]=tr.data
        else:
            if len(tr.data)<npts_avg or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_avg:
                continue
            data_avg[j,:]=tr.data

        if calc_type=="xcorr":
            stnm_avg[stnm_idx[tr.stats.station+"."+tr.stats.channel]]=True

        j+=1

    nsta_avg=j

def st2data_order():
    global nsta_avg

    i=0
    while i<len(st_avg):
        tr=st_avg[i]
        if bool_deci:
            if len(tr.data)<npts_avg_0 or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_avg_0:
                st_avg.pop(index=i)
                continue
        else:
            if len(tr.data)<npts_avg or type(tr.data)==ma.MaskedArray and tr.data.count()<npts_avg:
                st_avg.pop(index=i)
                continue

        stnm_avg[stnm_idx[tr.stats.station+"."+tr.stats.channel]]=True
        order[i]=stnm_idx[tr.stats.station+"."+tr.stats.channel]
        i+=1

    nsta_avg=len(st_avg)
    idx_order[:nsta_avg]=np.argsort(a=np.argsort(a=order[:nsta_avg]))
    for i in range(nsta_avg):
        tr=st_avg[i]
        if bool_deci:
            data_avg_0[idx_order[i],:]=tr.data
        else:
            data_avg[idx_order[i],:]=tr.data

def calc_day(idx_day):
    global st_avg,idx_triu_avg,idx_flat_avg

    for i in range(navg_day):
        if calc_type=="xcorr":
            stnm_avg[:]=0

        idx_avg=idx_day*navg_day+i
        start_avg=udt_cur+offset_avg+i*shift_avg
        st_avg=st_cur.slice(starttime=start_avg,endtime=start_avg+len_avg-delta)

        if st_step=="conti" and calc_type=="xcorr":
            st2data_order()
        else:
            st2data()

        idx_triu_avg=np.triu_indices(n=nsta_avg,k=1)
        if calc_type=="xcorr":
            idx_flat_avg=(stnm_avg[:,None] & stnm_avg)[idx_triu]

        pre_proc()
        calc_cov()
        if calc_type=="cohrnc":
            calc_cohrnc(idx_avg=idx_avg)
        elif calc_type=="xcorr":
            calc_xcorr(idx_avg=idx_avg)
        print_status(idx_avg=idx_avg)

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
    global st_cur,udt_cur,udt_next

    read_cur()
    for idx_day in range(nday):
        read_next()
        st_cur+=st_next.slice(endtime=udt_next+ovrlp_avg-delta)
        st_cur.merge(method=1)

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
