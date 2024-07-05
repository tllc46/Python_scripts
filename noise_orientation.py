import glob
from os.path import exists
from os import mkdir
from sys import stderr
from pprint import pprint

import numpy as np
import numpy.ma as ma

from obspy import read,read_inventory,UTCDateTime,Trace,Stream
from obspy.signal.util import next_pow_2
from obspy.geodetics import gps2dist_azimuth
from obspy.io.sac import SACTrace

sec_1d=86400

len_sub=3600
shift_sub=0.5*len_sub

n_sub_1d=int(sec_1d/shift_sub)

n_sub_tot=24
n_sub_shift=int(0.5*n_sub_tot)

sampling_rate_0=200
sampling_rate=10
factor=int(sampling_rate_0/sampling_rate)
delta=1/sampling_rate

npts_sub=sampling_rate*len_sub
n_freq=next_pow_2(i=npts_sub)
freq_min=0.015
freq_max=4.5
f_min_idx=int(np.ceil(freq_min*n_freq/sampling_rate))
f_max_idx=int(np.floor(freq_max*n_freq/sampling_rate))
n_mov_win=40
n_lag=2048
n_rfft=int(0.5*n_freq)+1

n_close=5

year=2020
start_julday=92
end_julday=96
start_1d=UTCDateTime(year=year,julday=start_julday)

net=read_inventory(path_or_file_object="stameta_NS",format="STATIONTXT",level="station")[0]
n_sta=len(net)

code_idx={}
pair_idx=np.empty(shape=(n_sta,n_close),dtype=int)
idx_pair=[]
sym_idx=np.zeros(shape=(n_sta,n_close,2),dtype=int)

avg_i_ls=[1]
n_sub_ls=[]
ZZ_ls=[]
NZ_ls=[]
EZ_ls=[]
n_ZZ_ls=[]
n_NZ_ls=[]
n_EZ_ls=[]

st_Z=Stream()
st_N=Stream()
st_E=Stream()

bool_sub_Z=np.full(shape=n_sta,fill_value=False)
bool_sub_N=np.full(shape=n_sta,fill_value=False)
bool_sub_E=np.full(shape=n_sta,fill_value=False)
bool_sub=np.full(shape=n_sta,fill_value=False)

spec_Z=np.empty(shape=(n_sta,n_rfft),dtype=complex)
spec_N=np.empty(shape=(n_sta,n_rfft),dtype=complex)
spec_E=np.empty(shape=(n_sta,n_rfft),dtype=complex)

sym_fill=np.full(shape=(n_sta,n_close),fill_value=False)

def code_idx_map():
    global code_idx,pair_idx,idx_pair,sym_idx
    sta_idx=np.zeros(shape=(n_sta,n_sta),dtype=int)

    for i in range(n_sta):
        code_idx[net[i].code]=i

    for i in range(n_sta):
        idx_pair.append([])
        i_idx=code_idx[net[i].code]
        gcarc_ls=[]
        for j in range(n_sta):
            if i!=j:
                gcarc_ls.append([gps2dist_azimuth(lat1=net[i].latitude,lon1=net[i].longitude,lat2=net[j].latitude,lon2=net[j].longitude)[0],net[j].code])
        gcarc_ls.sort()
        for j in range(n_close): #i번째 지진계와 가장 가까운 n_close개 지진계만
            j_idx=code_idx[gcarc_ls[j][1]]
            sta_idx[i_idx,j_idx]=j+1
            idx_pair[-1].append(gcarc_ls[j][1])

    for i in range(n_sta): #대칭 검사
        for j in range(i+1,n_sta):
            if sta_idx[i,j] and sta_idx[j,i]:
                sym_idx[i,sta_idx[i,j]-1]=[j+1,sta_idx[j,i]]

    for i in range(n_sta):
        for j in range(n_sta):
            if sta_idx[i,j]:
                pair_idx[i,sta_idx[i,j]-1]=j+1

def new_avg():
    global n_avg

    new_avg_i=avg_i_ls[-1]

    print(f"{new_avg_i:02} average window start",file=stderr)

    if not exists(path=f"/media/tllc46/data01/Orientation/{new_avg_i:02}"):
        mkdir(path=f"/media/tllc46/data01/Orientation/{new_avg_i:02}")
    for i in ["Z","N","E"]:
        if not exists(path=f"/media/tllc46/data01/Orientation/{new_avg_i:02}/{i}Z"):
            mkdir(path=f"/media/tllc46/data01/Orientation/{new_avg_i:02}/{i}Z")

    n_sub_ls.append(0)
    n_avg=len(n_sub_ls)
    avg_i_ls.append(new_avg_i+1)
    ZZ_ls.append(np.zeros(shape=(n_sta,n_close,2*n_lag+1)))
    NZ_ls.append(np.zeros(shape=(n_sta,n_close,2*n_lag+1)))
    EZ_ls.append(np.zeros(shape=(n_sta,n_close,2*n_lag+1)))
    n_ZZ_ls.append(np.zeros(shape=(n_sta,n_close),dtype=int))
    n_NZ_ls.append(np.zeros(shape=(n_sta,n_close),dtype=int))
    n_EZ_ls.append(np.zeros(shape=(n_sta,n_close),dtype=int))

def read_data(julday,code,component):
    pathname=f"/home/tllc46/NSNAS/{code}/HH{component}/*.{year}.{julday:03}"
    if not glob.glob(pathname=pathname):
        print(f"{code} {component} no data",file=stderr)
        return None
    else:
        st_temp=read(pathname_or_url=pathname,format="MSEED",header_byteorder=">")
        st_temp.merge(method=1)
        tr_temp=st_temp[0]
        if type(tr_temp.data)==ma.MaskedArray:
            print(f"{code} {component} discontinous data",file=stderr)
            return None
        if tr_temp.stats.component=="X":
            tr_temp.stats.component="E"
        elif tr_temp.stats.component=="Y":
            tr_temp.stats.component="N"
        tr_temp.decimate(factor=factor)
        tr_temp.filter(type="bandpass",freqmin=freq_min,freqmax=freq_max,zerophase=True)

        return tr_temp

def new_julday(julday):
    global st_Z,st_N,st_E

    print("reading data start",file=stderr)

    for i in range(n_sta):
        #Z
        tr_temp=read_data(julday=julday,code=net[i].code,component="Z")
        if not tr_temp:
            print(f"no {net[i].code} Z, skip reading N, E",file=stderr)
            continue
        else:
            st_Z+=tr_temp

        #N
        tr_temp=read_data(julday=julday,code=net[i].code,component="N")
        if not tr_temp:
            print(f"no {net[i].code} N, skip reading E",file=stderr)
            continue
        else:
            st_N+=tr_temp

        #E
        tr_temp=read_data(julday=julday,code=net[i].code,component="E")
        if tr_temp:
            st_E+=tr_temp

    print("reading data end\n",file=stderr)

    st_Z.merge(method=1)
    st_N.merge(method=1)
    st_E.merge(method=1)

def gap_check():
    global n_Z,n_N,n_E

    print("gap check",file=stderr)

    for st_i in [st_sub_Z,st_sub_N,st_sub_E]:
        j=0
        while j<len(st_i):
            tr_j=st_i[j]
            if type(tr_j.data)==ma.MaskedArray and tr_j.data.count()<0.5*len_sub*sampling_rate or len(tr_j.data)<0.5*len_sub*sampling_rate:
                print(f"{tr_j.stats.station} {tr_j.stats.component} gap is too large",file=stderr)
                st_i.pop(index=j)
                continue
            if type(tr_j.data)==ma.MaskedArray:
                tr_j[j].data.filled(fill_value=0)
            j+=1

    n_Z=len(st_sub_Z)
    n_N=len(st_sub_N)
    n_E=len(st_sub_E)

def sanity_check():
    global bool_sub,bool_sub_Z,bool_sub_N,bool_sub_E
    global spec_Z,spec_N,spec_E

    bool_sub_Z[:]=False
    bool_sub_N[:]=False
    bool_sub_E[:]=False

    #Z
    for i in range(n_Z):
        idx=code_idx[st_sub_Z[i].stats.station]
        bool_sub_Z[idx]=True

    #N
    for i in range(n_N):
        idx=code_idx[st_sub_N[i].stats.station]
        bool_sub_N[idx]=True

    #E
    for i in range(n_E):
        idx=code_idx[st_sub_E[i].stats.station]
        bool_sub_E[idx]=True

    bool_sub[:]=bool_sub_Z & bool_sub_N & bool_sub_E

def fourier_trans(data):
    spectrum=np.fft.rfft(a=np.sign(data),n=n_freq) #1비트 정규화, Fourier 변환
    norm=np.convolve(a=abs(spectrum),v=np.ones(shape=2*n_mov_win+1),mode="valid") #npwing 길이의 이동 평균
    result=np.zeros(shape=int(0.5*n_freq)+1,dtype=complex)
    result[f_min_idx:f_max_idx+1]=spectrum[f_min_idx:f_max_idx+1]/norm[f_min_idx-n_mov_win:f_max_idx+1-n_mov_win]*(2*n_mov_win+1) #spectral whitening
    return result

def spectra():
    global spec_Z,spec_N,spec_E

    #Z
    for i in range(n_Z):
        idx=code_idx[st_sub_Z[i].stats.station]
        if bool_sub_Z[idx]:
            spec_Z[idx]=fourier_trans(data=st_sub_Z[i].data)

    #N
    for i in range(n_N):
        idx=code_idx[st_sub_N[i].stats.station]
        if bool_sub[idx]:
            spec_N[idx]=fourier_trans(data=st_sub_N[i].data)

    #E
    for i in range(n_E):
        idx=code_idx[st_sub_E[i].stats.station]
        if bool_sub_E[idx]:
            spec_E[idx]=fourier_trans(data=st_sub_E[i].data)

def cross_correlate(spectrum1,spectrum2):
    norm=abs(spectrum1)*abs(spectrum2) #정규화 인자
    norm[0<norm]=1/norm[0<norm]
    corr_spec=spectrum1*np.conjugate(spectrum2)*norm
    corr=np.fft.irfft(a=corr_spec) #역 Fourier 변환
    corr=np.real(val=np.concatenate((corr[-n_lag:],corr[:n_lag+1]))) #time lag까지만 자르다
    return corr

def stacking():
    global ZZ_ls,NZ_ls,EZ_ls
    global n_ZZ_ls,n_NZ_ls,n_EZ_ls
    global sym_fill

    sym_fill[:,:]=False

    for i in range(n_sta): #행 방향
        for j in range(n_close): #열 방향
            idx=pair_idx[i,j]-1
            if bool_sub_Z[idx]:

                #Z
                if bool_sub[i]:
                    if not sym_fill[i,j]:
                        corr=cross_correlate(spectrum1=spec_Z[i],spectrum2=spec_Z[idx])
                        for k in range(n_avg):
                            ZZ_ls[k][i,j]+=corr
                            n_ZZ_ls[k][i,j]+=1

                    i_idx=sym_idx[i,j,0]-1
                    j_idx=sym_idx[i,j,1]-1
                    p_idx=pair_idx[i_idx,j_idx]-1
                    if sym_idx[i,j].any() and bool_sub[i_idx] and bool_sub_Z[p_idx] and not sym_fill[i_idx,j_idx]:
                        sym_fill[i_idx,j_idx]=True
                        for k in range(n_avg):
                            ZZ_ls[k][i_idx,j_idx]+=corr[::-1]
                            n_ZZ_ls[k][i_idx,j_idx]+=1

                #N
                if bool_sub[i]:
                    corr=cross_correlate(spectrum1=spec_N[i],spectrum2=spec_Z[idx])
                    for k in range(n_avg):
                        NZ_ls[k][i,j]+=corr
                        n_NZ_ls[k][i,j]+=1

                #E
                if bool_sub[i]:
                    corr=cross_correlate(spectrum1=spec_E[i],spectrum2=spec_Z[idx])
                    for k in range(n_avg):
                        EZ_ls[k][i,j]+=corr
                        n_EZ_ls[k][i,j]+=1

def old_jul_day():
    global start_1d

    print(f"{i:03} Julian day end",file=stderr)

    start_1d+=sec_1d
    st_Z.trim(starttime=start_1d)
    st_N.trim(starttime=start_1d)
    st_E.trim(starttime=start_1d)

def old_avg():
    global n_avg

    old_avg_i=avg_i_ls.pop(0)

    print(f"{old_avg_i:02} average window end\n",file=stderr)

    n_sub_ls.pop(0)
    n_avg=len(n_sub_ls)
    avg_ZZ=ZZ_ls.pop(0)
    avg_NZ=NZ_ls.pop(0)
    avg_EZ=EZ_ls.pop(0)
    n_sub_ZZ=n_ZZ_ls.pop(0)
    n_sub_NZ=n_NZ_ls.pop(0)
    n_sub_EZ=n_EZ_ls.pop(0)

    for i in range(n_sta):
        for j in range(n_close):
            #Z
            if n_sub_ZZ[i,j]:
                data=avg_ZZ[i,j]/n_sub_ZZ[i,j]
                sac=SACTrace.from_obspy_trace(trace=Trace(data=data,header={"delta":delta}))
                sac.user0=n_sub_ZZ[i,j]
                path=f"/media/tllc46/data01/Orientation/{old_avg_i:02}/ZZ/{net[i].code}-{idx_pair[i][j]}.SAC"
                sac.write(dest=path)

            #N
            if n_sub_NZ[i,j]:
                data=avg_NZ[i,j]/n_sub_NZ[i,j]
                sac=SACTrace.from_obspy_trace(trace=Trace(data=data,header={"delta":delta}))
                sac.user0=n_sub_NZ[i,j]
                path=f"/media/tllc46/data01/Orientation/{old_avg_i:02}/NZ/{net[i].code}-{idx_pair[i][j]}.SAC"
                sac.write(dest=path)

            #E
            if n_sub_EZ[i,j]:
                data=avg_EZ[i,j]/n_sub_EZ[i,j]
                sac=SACTrace.from_obspy_trace(trace=Trace(data=data,header={"delta":delta}))
                sac.user0=n_sub_EZ[i,j]
                path=f"/media/tllc46/data01/Orientation/{old_avg_i:02}/EZ/{net[i].code}-{idx_pair[i][j]}.SAC"
                sac.write(dest=path)

#main
code_idx_map()
new_julday(julday=start_julday)
new_avg()

st_Z.trim(starttime=start_1d,pad=True)
st_N.trim(starttime=start_1d,pad=True)
st_E.trim(starttime=start_1d,pad=True)

for i in range(start_julday,end_julday+1):
    print(f"{i:03} Julian day starting...",file=stderr)
    new_julday(julday=i+1)
    st_Z.trim(endtime=start_1d+2*sec_1d-delta,pad=True)

    for j in range(n_sub_1d):
        print(f"{j+1:02}/{n_sub_1d} sub window starting...",file=stderr)

        if n_sub_ls[-1]==n_sub_shift:
            new_avg()

        start_sub=start_1d+j*shift_sub
        st_sub_Z=st_Z.slice(starttime=start_sub,endtime=start_sub+len_sub-delta)
        st_sub_N=st_N.slice(starttime=start_sub,endtime=start_sub+len_sub-delta)
        st_sub_E=st_E.slice(starttime=start_sub,endtime=start_sub+len_sub-delta)

        gap_check()
        sanity_check()
        spectra()
        stacking()

        for k in range(n_avg):
            n_sub_ls[k]+=1

        if n_sub_ls[0]==n_sub_tot:
            old_avg()

        print(f"{j+1:02}/{n_sub_1d} sub window end\n",file=stderr)

    old_jul_day()
