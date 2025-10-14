import sys

from datetime import datetime
from importlib import import_module

import numpy as np
import pandas as pd

sys.path.append("/home/tllc46/Aso/tstep")
sys.path.append("/home/tllc46/Aso/xcorr/params")

mdl_t=import_module(name="t06")

def init(date_start,date_end):
    global dt_start,dt_end
    dt_start=datetime.strptime(date_start,"%Y-%m-%d")
    dt_end=datetime.strptime(date_end,"%Y-%m-%d")
    mdl_t.find_range(date_start=date_start,date_end=date_end)

def volc_eq(ax):
    df=pd.read_csv(filepath_or_buffer="/home/tllc46/48NAS1/tllc46/Aso/eq_data/num_earthquakes.csv",parse_dates=[0])
    df=df[(dt_start<=df["date"]) & (df["date"]<dt_end)]
    y=df[["A","BH","BP","B"]].to_numpy().T
    labels=["A","hi freq B","harmonic B","unid B"]
    handles=[]
    bottom=np.zeros(shape=len(df),dtype=int)

    for i in ["A","BH","BP","B"]:
        handles.append(ax.bar(x=df["date"],height=df[i],bottom=bottom))
        bottom+=df[i]

    ax.set_ylabel(ylabel="number of volcanic\nearthquakes (/day)")
    ax.set_yscale(value="log")
    ax.legend(handles=handles,labels=labels,loc="upper left")

def volc_tr(ax):
    df=pd.read_csv(filepath_or_buffer="/home/tllc46/48NAS1/tllc46/Aso/eq_data/num_earthquakes.csv",parse_dates=[0])
    df=df[(dt_start<=df["date"]) & (df["date"]<dt_end)]

    ax.bar(x=df["date"],height=df["T"],color="blue")

    ax.set_ylabel(ylabel="number of volcanic\ntremors (/day)")
    ax.set_yscale(value="log")

def conti_tr(ax):
    df=pd.read_csv(filepath_or_buffer="/home/tllc46/48NAS1/tllc46/Aso/eq_data/conti_tremor.csv",sep=" ",names=["date_start","date_end"],parse_dates=[0,1])
    df=df[(dt_start<df["date_end"]) & (df["date_start"]<dt_end)]
    df["width"]=df["date_end"]-df["date_start"]

    ax.broken_barh(xranges=df[["date_start","width"]].to_numpy(),yrange=(0,1),facecolors="blue")

    ax.set_yticks(ticks=[])

def iso_tr(ax):
    df=pd.read_csv(filepath_or_buffer="/home/tllc46/48NAS1/tllc46/Aso/eq_data/num_earthquakes.csv",parse_dates=[0])
    df=df[(dt_start<=df["date"]) & (df["date"]<dt_end)]

    ax.bar(x=df["date"],height=df["TK"],color="blue")

    ax.set_ylabel(ylabel="number of isolated\ntremors (/day)")
    ax.set_yscale(value="log")

def amp(ax,**kwargs):
    navg=len(mdl_t.ndt_range)
    amp=np.empty(shape=0)
    for i in mdl_t.names:
        amp_unit=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/spgm/V.ASOA.U/amp/"+i+".npz")
        amp_unit=amp_unit["amp"] #(navg,)
        amp=np.hstack(tup=(amp,amp_unit))
    amp=amp[mdl_t.idx_start:mdl_t.idx_start+navg]

    ax.plot(mdl_t.ndt_range,amp,color="black",**kwargs)

    ax.set_ylabel(ylabel="20m mean amplitude (nm/s)")
    ax.set_yscale(value="log")
    ax.set_ylim(bottom=10)

def spgm(ax,fig,**kwargs):
    mdl_x=import_module(name="x195")

    sampling_rate=100 #[Hz]
    delta=1/sampling_rate

    len_sub=mdl_x.len_sub
    npts_sub=len_sub*sampling_rate

    fmin_offset=0.1 #[Hz]
    fmin=1 #[Hz]
    fmax=10 #[Hz]
    idx_offset=int(np.floor(fmin_offset*len_sub))
    idx_fmin=int(np.floor(fmin*len_sub))
    idx_fmax=int(np.ceil(fmax*len_sub))+1
    nfreq=idx_fmax-idx_fmin
    frequency=np.fft.rfftfreq(n=npts_sub,d=delta)[idx_fmin:idx_fmax]
    idx_fmin-=idx_offset
    idx_fmax-=idx_offset

    navg=len(mdl_t.ndt_range)
    spgm=np.empty(shape=(nfreq,0))
    for i in mdl_t.names:
        spgm_unit=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/spgm/V.ASOA.U/foo1/"+i+".npz")
        spgm_unit=spgm_unit["spgm"][idx_fmin:idx_fmax] #(nfreq,navg)
        spgm=np.hstack(tup=(spgm,spgm_unit))
    spgm=spgm[:,mdl_t.idx_start:mdl_t.idx_start+navg]

    quadmesh=ax.pcolormesh(mdl_t.ndt_range,frequency,spgm,norm="log",cmap="rainbow",vmin=1e4,vmax=1e6)
    ax.axhline(y=3,color="black",linestyle="dashed")
    ax.axhline(y=6,color="black",linestyle="dashed")

    colorbar=fig.colorbar(mappable=quadmesh,ax=ax,**kwargs)
    colorbar.set_label(label="spectrum (nm/s)")

    ax.set_yscale(value="log")
    ax.set_ylabel(ylabel="frequency (Hz)")
    ax.yaxis.set_major_formatter(formatter="{x:g}")
