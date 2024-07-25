#obspy/signal/spectral_estimation.py 참고

import sys

import numpy as np
import numpy.ma as ma
from scipy.signal import ShortTimeFFT
import matplotlib.pyplot as plt

from obspy.signal.util import prev_pow_2
from obspy.signal.spectral_estimation import fft_taper

dtiny=np.finfo(dtype=float).tiny

sampling_rate=200
db_bins=(-200,-50,1)
ppsd_length=3600
overlap=0.5
period_smoothing_width_octaves=1
period_step_octaves=0.125

#13 segments overlapping by 75%
#1 segment + 25% * 12 segments
nfft=ppsd_length*sampling_rate*4/16
nfft=prev_pow_2(i=nfft)

nlap=int(0.75*nfft) #75% overlap
seg_len=int(sampling_rate*ppsd_length)

SFT=ShortTimeFFT(win=fft_taper(data=np.ones(shape=nfft)),hop=nfft-nlap,fs=sampling_rate,scale_to="psd",phase_shift=None)
_,p_lb=SFT.lower_border_end
_,p_ub=SFT.upper_border_begin(n=seg_len)
frequency_lim=(sampling_rate/nfft,0.5*sampling_rate)

num_frequency_bins=int((np.log2(nfft)-1)/period_step_octaves)+1
center_exp=np.arange(stop=num_frequency_bins)*period_step_octaves
edge_exp=center_exp+0.5*period_step_octaves
edge_exp=np.insert(arr=edge_exp,obj=0,values=-0.5*period_step_octaves)
frequency_bin_edges=frequency_lim[0]*np.power(2,edge_exp)

num_db_bins=int((db_bins[1]-db_bins[0])/db_bins[2])
db_bin_edges=np.linspace(start=db_bins[0],stop=db_bins[1],num=num_db_bins+1)

step=ppsd_length*(1-overlap)

stnm=sys.argv[1]

start_date="2013-10-01"
end_date="2015-10-31"

model_file="/home/tllc46/anaconda3/envs/seis/lib/python3.11/site-packages/obspy/signal/data/noise_models.npz"
NpzFile=np.load(file=model_file)
model_frequencies=1/NpzFile["model_periods"]
nhnm=NpzFile["high_noise"]
nlnm=NpzFile["low_noise"]

def calculate_psd():
    from sys import stderr
    from glob import glob

    from obspy import read,UTCDateTime
    from obspy.signal.invsim import evalresp

    sec_1d=86400
    udt_start_date=UTCDateTime(start_date)
    udt_end_date=UTCDateTime(end_date)
    start_1d=udt_start_date
    n_day=int((udt_end_date-udt_start_date)/sec_1d)+1
    n_seg_1d=int(sec_1d/step)
    n_seg=n_day*n_seg_1d

    left_exp=center_exp-0.5*period_smoothing_width_octaves
    right_exp=center_exp+0.5*period_smoothing_width_octaves
    psd_frequency_order=np.arange(start=1,stop=int(nfft/2)+1)
    psd_frequency_log2=np.log2(psd_frequency_order)
    left_idx=np.searchsorted(a=psd_frequency_log2,v=left_exp)
    right_idx=np.searchsorted(a=psd_frequency_log2,v=right_exp,side="right")

    delta=1/sampling_rate
    scaling_factor=2

    filename=f"/media/tllc46/data01/Jeju/BB_RESP/RESP.05.{stnm}..HHZ"
    resp=evalresp(t_samp=delta,nfft=nfft,filename=filename,date=UTCDateTime("2014-06-01"),units="ACC")
    resp=resp[1:]
    respamp=abs(resp*np.conj(resp))

    seg_idx=0
    psd_values=np.empty(shape=(n_seg,num_frequency_bins))

    def read_stream(date):
        pathname=f"/home/tllc46/48NAS2/symbolic.jeju/05.{stnm}/HHZ/05.{stnm}.HHZ.{date.year}.{date.julday:03}*"
        if glob(pathname=pathname):
            st=read(pathname_or_url=pathname,format="MSEED")
            l=0
            while l<len(st):
                if st[l].stats.sampling_rate!=200:
                    st.pop(index=l)
                    continue
                l+=1
            st.merge(method=1)
            return st
        else:
            return None

    def next_day():
        nonlocal st1,start_1d
        print(start_1d.strftime(format="%Y-%m-%d"),"| ending...",file=stderr)
        print("-------------------------------",file=stderr)
        st1=st2
        start_1d+=sec_1d

    #main loop
    st1=read_stream(date=start_1d)
    for i in range(n_day):
        print(start_1d.strftime(format="%Y-%m-%d"),"| starting...",file=stderr)
        st2=read_stream(date=start_1d+sec_1d)
        if not st1:
            print(start_1d.strftime(format="%Y-%m-%d"),"| no data",file=stderr)
            psd_values[seg_idx:seg_idx+n_seg_1d,:]=999999
            seg_idx+=n_seg_1d
            next_day()
            continue

        if st2:
            st1+=st2
        st1.merge(method=1)
        st1.trim(starttime=start_1d,endtime=start_1d+sec_1d+ppsd_length*overlap-delta)

        print(start_1d.strftime(format="%Y-%m-%d"),"| main processing...",file=stderr)
        for j in range(n_seg_1d):
            start_seg=start_1d+j*step
            st_seg=st1.slice(starttime=start_seg,endtime=start_seg+ppsd_length-delta)

            if not st_seg:
                print(start_1d.strftime(format="%Y-%m-%d"),f"| {j+1:02}/{n_seg_1d} segment | no data",file=stderr)
                psd_values[seg_idx,:]=999999
                seg_idx+=1
                continue

            tr=st_seg[0]

            if len(tr.data)<seg_len or type(tr.data)==ma.MaskedArray and tr.data.count()<seg_len:
                print(start_1d.strftime(format="%Y-%m-%d"),f"| {j+1:02}/{n_seg_1d} segment | gap exists",file=stderr)
                psd_values[seg_idx,:]=999999
                seg_idx+=1
                continue

            result=SFT.stft_detrend(x=tr.data,detr="linear",p0=p_lb,p1=p_ub)[1:]
            result=np.conj(result)*result
            result[:-1]*=scaling_factor
            spec=np.mean(a=result,axis=1).real
            spec/=respamp
            spec[spec<dtiny]=dtiny
            spec=10*np.log10(spec)
            for i in range(num_frequency_bins):
                psd_values[seg_idx,i]=np.mean(a=spec[left_idx[i]:right_idx[i]])
            seg_idx+=1

        next_day()

    np.save(file=f"ppsd_results/{stnm}.npy",arr=psd_values)

def plot_histogram():
    from obspy.imaging.cm import pqlx

    max_percentage=30
    cmap=pqlx

    psd_values=np.load(file=f"ppsd_results/{stnm}.npy")
    psd_values=ma.masked_equal(x=psd_values,value=999999)
    psd_values=ma.compress_rows(a=psd_values)

    current_histogram_count=len(psd_values)

    psd_values[psd_values<db_bins[0]]=db_bins[0]
    psd_values[db_bins[1]<psd_values]=db_bins[1]

    current_histogram=np.empty(shape=(num_db_bins,num_frequency_bins))

    for i in range(num_frequency_bins):
        current_histogram[:,i]=np.histogram(a=psd_values[:,i],bins=db_bin_edges)[0]*100/current_histogram_count

    meshgrid=np.meshgrid(frequency_bin_edges,db_bin_edges)

    fig=plt.figure(figsize=(25.6,12.67))
    ax=fig.subplots()
    quadmesh=ax.pcolormesh(meshgrid[0],meshgrid[1],current_histogram,cmap=cmap)
    ax.plot(model_frequencies,nhnm,color="0.4",linewidth=2)
    ax.plot(model_frequencies,nlnm,color="0.4",linewidth=2)
    colorbar=fig.colorbar(mappable=quadmesh,ax=ax)
    colorbar.mappable.set_clim(vmin=0,vmax=max_percentage)
    colorbar.set_label(label="[%]")
    quadmesh.set_clim(vmin=0,vmax=max_percentage)
    ax.grid(visible=True,which="both")
    ax.set_xlabel(xlabel="Frequency [Hz]")
    ax.set_xscale(value="log")
    ax.set_xlim(left=frequency_lim[0],right=frequency_lim[1])
    ax.set_ylim(bottom=db_bins[0],top=db_bins[1])
    ax.set_ylabel(ylabel="Amplitude [$(m/s^2)^2/Hz$] [dB]")
    ax.xaxis.set_major_formatter(formatter="{x:g}")
    ax.set_title(label=stnm)
    fig.savefig(fname=f"ppsd_plots/{stnm}.png")

def plot_statistics():
    percentiles=[0,25,50,75,100]
    frequency_bin_centers=frequency_lim[0]*np.power(2,center_exp)
    db_bin_centers=0.5*(db_bin_edges[:-1]+db_bin_edges[1:])

    psd_values=np.load(file=f"ppsd_results/{stnm}.npy")
    psd_values=ma.masked_equal(x=psd_values,value=999999)
    psd_values=ma.compress_rows(a=psd_values)

    current_histogram=np.empty(shape=(num_db_bins,num_frequency_bins),dtype=int)

    for i in range(num_frequency_bins):
        current_histogram[:,i]=np.histogram(a=psd_values[:,i],bins=db_bin_edges)[0]

    percentile_values=np.percentile(a=psd_values,q=percentiles,axis=0)
    mode=db_bin_centers[np.argmax(a=current_histogram,axis=0)]
    mean=np.mean(a=psd_values,axis=0)

    fig=plt.figure(figsize=(25.6,12.67))
    ax=fig.subplots()
    for i in range(len(percentiles)-1):
        ax.plot(frequency_bin_centers,percentile_values[i],color="green")
    line1=ax.plot(frequency_bin_centers,percentile_values[i+1],color="green")[0]
    line2=ax.plot(frequency_bin_centers,mode,color="blue")[0]
    line3=ax.plot(frequency_bin_centers,mean,color="red")[0]
    ax.plot(model_frequencies,nhnm,color="0.4",linewidth=2)
    ax.plot(model_frequencies,nlnm,color="0.4",linewidth=2)
    ax.legend(handles=[line1,line2,line3],labels=["percentiles","mode","mean"])
    ax.grid(visible=True,which="both")
    ax.set_xlabel(xlabel="Frequency [Hz]")
    ax.set_xscale(value="log")
    ax.set_xlim(left=frequency_lim[0],right=frequency_lim[1])
    ax.set_ylim(bottom=db_bins[0],top=db_bins[1])
    ax.set_ylabel(ylabel="Amplitude [$(m/s^2)^2/Hz$] [dB]")
    ax.xaxis.set_major_formatter(formatter="{x:g}")
    ax.set_title(label=stnm)
    fig.savefig(fname=f"ppsd_plots/{stnm}.png")

def plot_spectrogram():
    np_start_date=np.datetime64(start_date)
    np_end_date=np.datetime64(end_date)
    np_step=np.timedelta64(int(step),"s")
    xedges=np.arange(start=np_start_date,stop=np_end_date+np.timedelta64(1,"D")+np_step,step=np_step)
    meshgrid_x,meshgrid_y=np.meshgrid(xedges,frequency_bin_edges)

    psd_values=np.load(file=f"ppsd_results/{stnm}.npy")
    psd_values=ma.masked_equal(x=psd_values,value=999999)
    psd_values=psd_values.T

    fig=plt.figure(figsize=(25.6,12.67))
    ax=fig.subplots()
    quadmesh=ax.pcolormesh(meshgrid_x,meshgrid_y,psd_values,cmap="rainbow")
    colorbar=fig.colorbar(mappable=quadmesh,ax=ax)
    colorbar.set_label(label="Amplitude [$(m/s^2)^2/Hz$] [dB]")
    ax.set_ylabel(ylabel="Frequency [Hz]")
    ax.set_yscale(value="log")
    ax.set_xlim(left=xedges[0],right=xedges[-1])
    ax.set_ylim(bottom=frequency_lim[0],top=frequency_lim[1])
    fig.savefig(fname=f"ppsd_plots/{stnm}.png")

#main
calculate_psd()
