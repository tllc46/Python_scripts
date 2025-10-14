from datetime import datetime

import pandas as pd

df=pd.read_csv(filepath_or_buffer="/home/tllc46/48NAS1/tllc46/Aso/vo_data/6_yudamari_info.csv",sep=" ",names=["date","area","T_surf","T_wall","T_floor","under","level"],parse_dates=[0])

def init(date_start,date_end):
    global df

    dt_start=datetime.strptime(date_start,"%Y-%m-%d")
    dt_end=datetime.strptime(date_end,"%Y-%m-%d")
    df=df[(dt_start<=df["date"]) & (df["date"]<dt_end)]

def area(ax):
    ax.scatter(x=df["date"],y=df["area"],color="black")
    ax.set_ylabel(ylabel="caldera area\n(proportion)")

def level(ax):
    ax.scatter(x=df["date"],y=df["level"],color="blue")
    ax.set_ylabel(ylabel="caldera\nwater level (m)")

def surf_temp(ax):
    ax.scatter(x=df["date"],y=df["T_surf"],color="blue")
    ax.set_ylabel(ylabel="caldera surface\ntemperature (°C)")

def fum_temp(ax):
    pcollect_1=ax.scatter(x=df["date"],y=df["T_wall"],color="green",marker="s")
    pcollect_2=ax.scatter(x=df.loc[df["under"]==0,"date"],y=df.loc[df["under"]==0,"T_floor"],color="blue")
    pcollect_3=ax.scatter(x=df.loc[df["under"]==1,"date"],y=df.loc[df["under"]==1,"T_floor"],color="red",marker="v")
    ax.set_ylabel(ylabel="fumarole\ntemperature (°C)")
    ax.legend(handles=[pcollect_1,pcollect_2,pcollect_3],labels=["wall","floor (explicit)","floor (up lim)"],loc="upper left")
