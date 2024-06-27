from datetime import datetime,timedelta
from os.path import exists
from os import mkdir

import numpy as np
import numpy.ma as ma
from scipy.signal import correlate

from obspy import read,read_inventory,UTCDateTime,Trace,Stream
from obspy.signal.util import next_pow_2
from obspy.geodetics import gps2dist_azimuth
from obspy.io.sac import SACTrace

subwin_start=UTCDateTime("2020-04-01")
endtime=UTCDateTime("2020-04-05")
subwin_length=3600
subwin_overlap=0.5
sampling_rate=10
delta=1/sampling_rate
subwin_end=subwin_start+subwin_length
nfreq=next_pow_2(i=subwin_length*sampling_rate)
nlag=2048
npwing=40
freqmin=0.015
freqmax=4.5
nmin=int(np.ceil(freqmin*nfreq/sampling_rate))
nmax=int(np.floor(freqmax*nfreq/sampling_rate))
win_length=3600*24
win_overlap=0.5
win_no=0
max_subwin=int(win_length/(subwin_length*(1-subwin_overlap)))
next_subwin=int(win_length*(1-subwin_overlap)/(subwin_length*(1-subwin_overlap)))

julday_ls=0
nsubwin_ls=[]
C_ZZ={}
C_NZ={}
C_EZ={}
st=Stream()

net=read_inventory(path_or_file_object="stameta_NS",format="STATIONTXT",level="station")[0] #network
code_idx={}
for i in range(len(net)):
    code_idx[net[i].code]=i

r_s_idx=np.zeros(shape=(len(net),len(net)),dtype=int)

for i in range(len(net)):
    r_idx=code_idx[net[i].code] #지진계
    gcarc_ls=[]
    for j in range(len(net)):
        if i!=j:
            gcarc_ls.append([gps2dist_azimuth(lat1=net[i].latitude,lon1=net[i].longitude,lat2=net[j].latitude,lon2=net[j].longitude)[0],net[j].code])
    gcarc_ls.sort()
    for j in range(20): #i번째 지진계와 가장 가까운 20개 지진계만
        s_idx=code_idx[gcarc_ls[j][1]] #source
        r_s_idx[r_idx,s_idx]=j+1

def fourier_trans(data):
    if type(data)==ma.MaskedArray:
        data=data.filled(fill_value=0) #공백은 0으로 채운다
    spectrum=np.fft.rfft(a=np.sign(data),n=nfreq) #1비트 정규화, Fourier 변환
    norm=correlate(in1=abs(spectrum),in2=np.ones(shape=2*npwing+1),mode="valid") #npwing 길이의 이동 평균
    result=np.zeros(shape=int(nfreq/2+1),dtype=complex)
    result[nmin:nmax+1]=spectrum[nmin:nmax+1]/norm[nmin-npwing:nmax+1-npwing]*(2*npwing+1) #spectral whitening
    return result

def cross_correlate(spectrum1,spectrum2):
    norm=np.abs(spectrum1)*np.abs(spectrum2) #정규화 인자
    norm[0<norm]=1/norm[0<norm]
    corr_spec=spectrum1*np.conjugate(spectrum2)*norm
    corr=np.fft.irfft(a=corr_spec) #역 Fourier 변환
    corr=np.real(np.concatenate((corr[-nlag:],corr[:nlag+1]))) #time lag까지만 자르다
    return corr

def preprocess(code,component):
    if not exists(path=f"/home/tllc46/NSNAS/{net[i].code}/HH{component}/NS.{net[i].code}.HH{component}.{subwin_end.year}.{subwin_end.julday:03}"): #miniseed 파일이 없다
        return None #제외
    else: #miniseed 파일 존재
        st_temp=read(pathname_or_url=f"/home/tllc46/NSNAS/{code}/HH{component}/*{subwin_end.year}.{subwin_end.julday:03}",format="MSEED",header_byteorder=">")
        print(f"done reading {code}")
        st_temp.merge(method=1)
        tr_temp=st_temp[0]
        if type(tr_temp.data)==ma.MaskedArray: #공백 존재
            return None #제외
        tr_temp.decimate(factor=int(tr_temp.stats.sampling_rate/sampling_rate),no_filter=True) #decimate
        tr_temp.filter(type="bandpass",freqmin=freqmin,freqmax=freqmax,zerophase=True) #대역 통과 필터
        if tr_temp.stats.component=="X": #성분 이름 수정
            tr_temp.stats.component="E"
        elif tr_temp.stats.component=="Y":
            tr_temp.stats.component="N"
        return tr_temp

def new_julday():
    global st_Z,st_N,st_E
    for i in range(len(net)): #각 지진계마다
        tr_temp=preprocess(code=net[i].code,component="Z")
        if not tr_temp: #Z축 성분이 존재하지 않다
            continue
        else:
            st_Z+=tr_temp

        tr_temp=preprocess(code=net[i].code,componenet="N")
        if not tr_temp: #N축 성분이 존재하지 않다
            continue
        else:
            st_N+=tr_temp

        tr_temp=preprocess(code=net[i].code,component="E")
        if tr_temp:
            st_E+=tr_temp
            
    st_Z.merge()
    st_Z.trim(starttime=subwin_start,endtime=UTCDateTime(year=subwin_end.year,julday=subwin_end.julday)+86400-delta,pad=True)
    st_N.merge()
    st_N.trim(starttime=subwin_start,endtime=UTCDateTime(year=subwin_end.year,julday=subwin_end.julday)+86400-delta,pad=True)
    st_E.merge()
    st_E.trim(starttime=subwin_start,endtime=UTCDateTime(year=subwin_end.year,julday=subwin_end.julday)+86400-delta,pad=True)

def C_new_win():
    global nsubwin_ls,C_ZZ,C_NZ,C_EZ
    if not exists(path=f"/media/tllc46/data01/Orientation/{win_no+1:02}"):
        mkdir(path=f"/media/tllc46/data01/Orientation/{win_no+1:02}")
    for i in ["Z","N","E"]:
        if not exists(path=f"/media/tllc46/data01/Orientation/{win_no+1:02}/{i}Z"):
            mkdir(path=f"/media/tllc46/data01/Orientation/{win_no+1:02}/{i}Z")
    nsubwin_ls.append(0)
    for i in C_ZZ:
        C_ZZ[i].append([np.zeros(shape=2*nlag+1),0])
    for i in C_NZ:
        C_NZ[i].append([np.zeros(shape=2*nlag+1),0])
        C_EZ[i].append([np.zeros(shape=2*nlag+1),0])

def C_old_win():
    global win_no,nsubwin_ls,C_ZZ,C_NZ,C_EZ
    win_no+=1
    nsubwin_ls.pop(0)
    for i in C_ZZ:
        subwin=C_ZZ[i].pop(0)
        if subwin[1]:
            data=subwin[0]/subwin[1]
            tr=Trace(data=data,header={"delta":delta})
            sac=SACTrace.from_obspy_trace(trace=tr)
            sac.user0=subwin[1]
            sac.write(dest=f"/media/tllc46/data01/Orientation/{win_no:02}/ZZ/{i}.SAC")
    for i in C_NZ:
        subwin=C_NZ[i].pop(0)
        if subwin[1]:
            data=subwin[0]/subwin[1]
            tr=Trace(data=data,header={"delta":delta})
            sac=SACTrace.from_obspy_trace(trace=tr)
            sac.user0=subwin[1]
            sac.write(dest=f"/media/tllc46/data01/Orientation/{win_no:02}/NZ/{i}.SAC")
        subwin=C_EZ[i].pop(0)
        if subwin[1]:
            data=subwin[0]/subwin[1]
            tr=Trace(data=data,header={"delta":delta})
            sac=SACTrace.from_obspy_trace(trace=tr)
            sac.user0=subwin[1]
            sac.write(dest=f"/media/tllc46/data01/Orientation/{win_no:02}/EZ/{i}.SAC")

C_new_win()

while subwin_end<=endtime:
    print(f"sub-window starts: {subwin_start} ~ {subwin_end}")
    #새 trace 추가
    if subwin_end.julday!=julday_ls: #새 Julian일 시작
        julday_ls=subwin_end.julday
        new_julday()

    #부구간으로만 자르다
    st_subwin_Z=st_Z.slice(starttime=subwin_start,endtime=subwin_end)
    st_subwin_N=st_N.slice(starttime=subwin_start,endtime=subwin_end)
    st_subwin_E=st_E.slice(starttime=subwin_start,endtime=subwin_end)

    #공백 검사
    for st_i in [st_subwin_Z,st_subwin_N,st_subwin_E]:
        j=0
        while j<len(st_i): #for all traces
            if type(st_i[j].data)==ma.MaskedArray and st_i[j].data.count()<subwin_length*sampling_rate*0.5: #공백이 부구간 길이의 50%보다 큰 경우
                st_i.pop(index=j)
                continue
            j+=1

    #sanity check
    for i in range(len(st_subwin_Z)):
        r_idx=code_idx[st_subwin_Z[i].stats.station]
        for j in range(i+1,len(st_subwin_Z)):
            s_idx=code_idx[st_subwin_Z[j].stats.station]
            if r_s_idx[r_idx,s_idx]:
                subwin_idx_Z[r_idx,r_s_idx[r_idx,s_idx]-1]=[i,j]
                sym+=1
            if r_s_idx[s_idx,r_idx]:
                subwin_idx_Z[s_idx,r_s_idx[s_idx,r_idx]-1]=[j,i]
                sym+=1
            if sym==2:
                sym_flag[s_idx,r_s_idx[s_idx,r_idx]-1]=True

    for i in range(len(st_subwin_N)):
        r_idx=code_idx[st_subwin_N[i].stats.station]
        for j in range(i+1,len(st_subwin_N)):
            s_idx=code_idx[st_subwin_N[j].stats.station]
            if r_s_idx[r_idx,s_idx]:
                subwin_idx_Z[r_idx,r_s_idx[r_idx,s_idx]-1]=[i,j]
            if r_s_idx[s_idx,r_idx]:
                subwin_idx_Z[s_idx,r_s_idx[s_idx,r_idx]-1]=[j,i]

                
    #cross correlate
    for i in range(len(vertical)):
        for j in range(i+1,len(vertical)):
            hor_exists=False
            st_i=st_subwin.select(station=vertical[i]) #E, N, Z order
            st_j=st_subwin.select(station=vertical[j])
            if len(st_i)==3 and f"{vertical[i]}-{vertical[j]}" in C_EZ: #cross correlate Ni-Zj, Ei-Zj
                hor_exists=True
                spec_Zj=fourier_trans(st_j[-1].data)
                spec_Ei=fourier_trans(st_i[0].data)
                corr=cross_correlate(spec_Ei,spec_Zj)
                for k in range(len(C_EZ[f"{vertical[i]}-{vertical[j]}"])):
                    C_EZ[f"{vertical[i]}-{vertical[j]}"][k][0]+=corr
                    C_EZ[f"{vertical[i]}-{vertical[j]}"][k][1]+=1
                spec_Ni=fourier_trans(st_i[1].data)
                corr=cross_correlate(spec_Ni,spec_Zj)
                for k in range(len(C_EZ[f"{vertical[i]}-{vertical[j]}"])):
                    C_NZ[f"{vertical[i]}-{vertical[j]}"][k][0]+=corr
                    C_NZ[f"{vertical[i]}-{vertical[j]}"][k][1]+=1

            if len(st_j)==3 and f"{vertical[j]}-{vertical[i]}" in C_EZ: #cross correlate Nj-Zi, Ej-Zi
                hor_exists=True
                spec_Zi=fourier_trans(st_i[-1].data)
                spec_Ej=fourier_trans(st_j[0].data)
                corr=cross_correlate(spec_Ej,spec_Zi)
                for k in range(len(C_EZ[f"{vertical[j]}-{vertical[i]}"])):
                    C_EZ[f"{vertical[j]}-{vertical[i]}"][k][0]+=corr
                    C_EZ[f"{vertical[j]}-{vertical[i]}"][k][1]+=1
                spec_Nj=fourier_trans(st_j[1].data)
                corr=cross_correlate(spec_Nj,spec_Zi)
                for k in range(len(C_EZ[f"{vertical[j]}-{vertical[i]}"])):
                    C_NZ[f"{vertical[j]}-{vertical[i]}"][k][0]+=corr
                    C_NZ[f"{vertical[j]}-{vertical[i]}"][k][1]+=1

            if hor_exists: #cross correlate Zi-Zj
                corr=cross_correlate(spec_Zi,spec_Zj)
                for k in range(len(C_ZZ[f"{vertical[i]}-{vertical[j]}"])):
                    C_ZZ[f"{vertical[i]}-{vertical[j]}"][k][0]+=corr
                    C_ZZ[f"{vertical[i]}-{vertical[j]}"][k][1]+=1
    
    print(f"sub-window ends: {subwin_start} ~ {subwin_end}")

    for i in range(len(nsubwin_ls)):
        nsubwin_ls[i]+=1

    if nsubwin_ls[0]==max_subwin:
        C_old_win()

    if nsubwin_ls[-1]==next_subwin:
        C_new_win()

    subwin_start+=subwin_length*(1-win_overlap)
    subwin_end+=subwin_length*(1-win_overlap)
