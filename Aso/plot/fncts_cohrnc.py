import sys
from importlib import import_module

import numpy as np
from matplotlib.axes import Axes

sys.path.append("/home/tllc46/Aso/tstep")
sys.path.append("/home/tllc46/Aso/xcorr/params")

mdl_t_day=import_module(name="t04")
mdl_t_avg=import_module(name="t06")
mdl_x=import_module(name="x195")

#sampling rate
sampling_rate=100 #[Hz]
delta=1/sampling_rate

#sub window
len_sub=mdl_x.len_sub
npts_sub=len_sub*sampling_rate

def init(date_start,date_end):
    mdl_t_day.find_range(date_start=date_start,date_end=date_end)
    mdl_t_avg.find_range(date_start=date_start,date_end=date_end)

def cohrnc(fig,ax,**kwargs):
    fmin_offset=0.1 #[Hz]
    fmin=1 #[Hz]
    fmax=10 #[Hz]
    idx_offset=int(np.floor(fmin_offset*len_sub))
    idx_fmin=int(np.floor(fmin*len_sub))
    idx_fmax=int(np.ceil(fmax*len_sub))+1
    nfreq=idx_fmax-idx_fmin
    frequency=np.fft.rfftfreq(n=npts_sub,d=delta)[idx_fmin:idx_fmax]
    frequency=np.append(arr=frequency,values=frequency[-1]+1/len_sub)-1/(2*len_sub)
    idx_fmin-=idx_offset
    idx_fmax-=idx_offset

    navg=len(mdl_t_day.ndt_range)

    cohrnc=np.empty(shape=(nfreq,0))
    for i in mdl_t_day.names:
        cohrnc_unit=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/t04/cohrnc/c197."+i+".npz")
        cohrnc_unit=cohrnc_unit["cohrnc"][idx_fmin:idx_fmax] #(nfreq,navg)
        cohrnc=np.hstack(tup=(cohrnc,cohrnc_unit))
    cohrnc=cohrnc[:,mdl_t_day.idx_start:mdl_t_day.idx_start+navg]

    ndt_range=np.append(arr=mdl_t_day.ndt_range,values=mdl_t_day.ndt_range[-1]+np.timedelta64(1,"D"))
    quadmesh=ax.pcolormesh(ndt_range,frequency,cohrnc,cmap="rainbow")
    ax.axhline(y=3,color="black",linestyle="dashed")
    ax.axhline(y=6,color="black",linestyle="dashed")

    colorbar=fig.colorbar(mappable=quadmesh,ax=ax,**kwargs)
    colorbar.set_label(label="coherence")

    ax.set_yscale(value="log")
    ax.set_ylabel(ylabel="frequency (Hz)")
    ax.yaxis.set_major_formatter(formatter="{x:g}")

def cohrnc_mean(ax,**kwargs):
    fmin_offset=0.1 #[Hz]
    fmin=3 #[Hz]
    fmax=6 #[Hz]
    idx_offset=int(np.floor(fmin_offset*len_sub))
    idx_fmin=int(np.floor(fmin*len_sub))
    idx_fmax=int(np.ceil(fmax*len_sub))+1
    nfreq=idx_fmax-idx_fmin
    frequency=np.fft.rfftfreq(n=npts_sub,d=delta)[idx_fmin:idx_fmax]
    idx_fmin-=idx_offset
    idx_fmax-=idx_offset

    navg=len(mdl_t_avg.ndt_range)

    cohrnc=np.empty(shape=(nfreq,0))
    for i in mdl_t_avg.names:
        cohrnc_unit=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/t06/cohrnc/c197."+i+".npz")
        cohrnc_unit=cohrnc_unit["cohrnc"][idx_fmin:idx_fmax] #(nfreq,navg)
        cohrnc=np.hstack(tup=(cohrnc,cohrnc_unit))
    cohrnc=cohrnc[:,mdl_t_avg.idx_start:mdl_t_avg.idx_start+navg]
    cohrnc_mean=np.mean(a=cohrnc,axis=0)

    ax.plot(mdl_t_avg.ndt_range,cohrnc_mean,**kwargs)

    ax.set_ylabel(ylabel="mean coherence\n(3-6Hz)")

def nsta_avg(ax):
    navg=len(mdl_t_avg.ndt_range)

    nsta_avg=np.empty(shape=0)
    for i in mdl_t_avg.names:
        nsta_unit=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/t06/cohrnc/c195."+i+".npz")
        nsta_unit=nsta_unit["nsta_avail"] #(navg,)
        nsta_avg=np.hstack(tup=(nsta_avg,nsta_unit))
    nsta_avg=nsta_avg[mdl_t_avg.idx_start:mdl_t_avg.idx_start+navg]

    ax.step(x=mdl_t_avg.ndt_range,y=nsta_avg,where="mid",color="blue")

    ax.set_ylabel(ylabel="no. of stations")
