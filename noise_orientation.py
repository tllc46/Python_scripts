from datetime import datetime,timedelta
from os.path import exists
from os import mkdir
from pprint import pprint

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
delta=0.1
subwin_end=subwin_start+subwin_length-delta
nfreq=next_pow_2(subwin_length*sampling_rate)
nlag=2048
npwing=40
freqmin=0.015
freqmax=4.5
nmin=int(np.ceil(freqmin*nfreq/sampling_rate))
nmax=int(np.floor(freqmax*nfreq/sampling_rate))
win_length=3600*24 #86400*30
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

for i in range(len(net)):
    gcarc_ls=[]
    for j in range(len(net)):
        if i!=j:
            gcarc_ls.append([gps2dist_azimuth(lat1=net[i].latitude,lon1=net[i].longitude,lat2=net[j].latitude,lon2=net[j].longitude)[0],net[j].code])
    gcarc_ls.sort()
    for j in range(20):
        C_NZ[f"{net[i].code}-{gcarc_ls[j][1]}"]=[]
        C_EZ[f"{net[i].code}-{gcarc_ls[j][1]}"]=[]
        if net[i].code<gcarc_ls[j][1] and f"{net[i].code}-{gcarc_ls[j][1]}" not in C_ZZ:
            C_ZZ[f"{net[i].code}-{gcarc_ls[j][1]}"]=[]
        elif gcarc_ls[j][1]<net[i].code and f"{gcarc_ls[j][1]}-{net[i].code}" not in C_ZZ:
            C_ZZ[f"{gcarc_ls[j][1]}-{net[i].code}"]=[]

def fourier_trans(data):
    if type(data)==ma.MaskedArray:
        data=data.filled(fill_value=0) #fill gaps with 0
    spectrum=np.fft.rfft(a=np.sign(data),n=nfreq) #1 bit normalize and fft
    norm=correlate(in1=np.absolute(spectrum),in2=np.ones(shape=2*npwing+1),mode="valid") #npwing moving window
    result=np.zeros(shape=int(nfreq/2+1),dtype=complex)
    result[nmin:nmax+1]=spectrum[nmin:nmax+1]/norm[nmin-npwing:nmax+1-npwing]*(2*npwing+1) #spectral whitening
    return result

def cross_correlate(spectrum1,spectrum2):
    norm=np.abs(spectrum1)*np.abs(spectrum2) #normalize factor
    norm[0<norm]=1/norm[0<norm]
    corr_spec=spectrum1*np.conjugate(spectrum2)*norm
    corr=np.fft.irfft(a=corr_spec) #ifft
    corr=np.real(np.concatenate((corr[-nlag:],corr[:nlag+1]))) #slice for time lag
    return corr

def new_julday():
    global julday_ls,st
    julday_ls=subwin_end.julday
    for i in range(len(net)): #for all stations
        for component in ["Z","N","E"]:
            if exists(path=f"/home/tllc46/NSNAS/{net[i].code}/HH{component}/NS.{net[i].code}.HH{component}.{subwin_end.year}.{subwin_end.julday:03}"): #miniseed file exists
                st_temp=read(pathname_or_url=f"/home/tllc46/NSNAS/{net[i].code}/HH{component}/*{subwin_end.year}.{subwin_end.julday:03}",format="MSEED",header_byteorder=">") #temporary stream for different sampling rate
                print(f"done reading {net[i].code}")
                st_temp.merge(method=1)
                tr_temp=st_temp[0]
                if type(tr_temp.data)==ma.MaskedArray: #inter-gap exists
                    continue #discard
                tr_temp.decimate(factor=int(tr_temp.stats.sampling_rate/sampling_rate),no_filter=True) #decimate
                tr_temp.filter(type="bandpass",freqmin=freqmin,freqmax=freqmax,zerophase=True) #filter
                if tr_temp.stats.component=="X": #correct for components
                    tr_temp.stats.component="E"
                elif tr_temp.stats.component=="Y":
                    tr_temp.stats.component="N"
                st+=tr_temp #from now, add this new stream to main stream
            else: #neither one of component exists, discard this station
                break
    st.merge() #merge again with old traces
    st.trim(starttime=subwin_start,endtime=UTCDateTime(datetime.strptime(f"{subwin_end.year} {subwin_end.julday}","%Y %j")+timedelta(days=1))-delta,pad=True) #append masked array if less than end of day

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
    vertical=[]
    print(f"sub-window starts: {subwin_start} ~ {subwin_end}")
    #add new trace
    if subwin_end.julday!=julday_ls: #new julian day starts
        new_julday()

    #slice for time window
    st_subwin=st.slice(starttime=subwin_start,endtime=subwin_end)

    #gap check
    i=0
    while i<len(st_subwin): #for all traces
        if type(st_subwin[i].data)==ma.MaskedArray and st_subwin[i].data.count()<subwin_length*sampling_rate*0.5: #gap exists and is larger than half of length of window
            st_subwin.pop(index=i)
        elif type(st_subwin[i].data)==ma.MaskedArray and subwin_length*sampling_rate*0.5<=st_subwin[i].data.count() or type(st_subwin[i].data)!=ma.MaskedArray: #gap doesn't exist or gap exists and is smaller than half of length of window 
            if st_subwin[i].stats.component=="Z" and st_subwin[i].stats.station not in vertical:
                vertical.append(st_subwin[i].stats.station)
            i+=1

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
