import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

path_save="/home/tllc46/48NAS1/tllc46/Aso/loc/tt/vel.npz"

topo=np.load(file="/home/tllc46/Downloads/topo/Aso_topo_10m.npz")
topo=topo["topo"] #(6000,6750)=(lat_topo,lon_topo) [m]
lat_topo=np.arange(stop=6000)/9000+(32.5835-1/9000)
lon_topo=np.arange(stop=6750)/9000+(130.7505-4/9000)
interp_topo=RegularGridInterpolator(points=(lat_topo,lon_topo),values=topo)

lon_orig=np.array(object=[130.95,131,131.05,131.1,131.15])
lat_orig=np.array(object=[32.8,32.85,32.9,32.95,33])
ele_orig=np.array(object=[-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0])
vel_orig=np.empty(shape=(5,5,13)) #(lon_orig,lat_orig,ele_orig)

for i in range(5): #lon_orig
    for j in range(5): #lat_orig
        filepath="/home/tllc46/Documents/Aso/loc/velocity/"+f"{lat_orig[j]:.2f}-{lon_orig[i]:.2f}"
        df=pd.read_csv(filepath_or_buffer=filepath,sep="\\s+",names=["depth","velocity"])
        vel_orig[i,j]=np.flip(m=df["velocity"].to_numpy())
interp_vel=RegularGridInterpolator(points=(lon_orig,lat_orig,ele_orig),values=vel_orig)

lon_min=130.95
lon_max=131.15
lon_res_inv=500
lon_res=1/lon_res_inv
lon_min_int=round(lon_min*lon_res_inv)
lon_max_int=round(lon_max*lon_res_inv)
lon=lon_res*np.arange(start=lon_min_int,stop=lon_max_int+1)
lon_num=lon_max_int+1-lon_min_int

lat_min=32.8
lat_max=33
lat_res_inv=500
lat_res=1/lat_res_inv
lat_min_int=round(lat_min*lat_res_inv)
lat_max_int=round(lat_max*lat_res_inv)
lat=lat_res*np.arange(start=lat_min_int,stop=lat_max_int+1)
lat_num=lat_max_int+1-lat_min_int

ele_min=-2.6 #[km]
ele_max=1.6 #[km]
ele_res_inv=25
ele_res=1/ele_res_inv
ele_min_int=round(ele_min*ele_res_inv)
ele_max_int=round(ele_max*ele_res_inv)
ele=ele_res*np.arange(start=ele_min_int,stop=ele_max_int+1)
ele_num=ele_max_int+1-ele_min_int

print(lon_num,lat_num,ele_num)
print(lon_res,lat_res,ele_res)

vel=np.zeros(shape=(lon_num,lat_num,ele_num)) #(lon,lat,ele)
idx=np.zeros(shape=(lon_num,lat_num,ele_num),dtype=bool) #(lon,lat,ele)

for i in range(lon_num): #lon
    for j in range(lat_num): #lat
        topo_ij=interp_topo(xi=(lat[j],lon[i]))
        topo_ij=ele_res*round(topo_ij*ele_res_inv/1000) #round to ele_res [km]
        for k in range(ele_num): #ele
            ele_k=ele[k] #[km]
            #print(lon[i],lat[j],ele_k-topo_ij)
            if topo_ij-6<=ele_k and ele_k<=topo_ij:
                vel[i,j,k]=interp_vel(xi=(lon[i],lat[j],ele_k-topo_ij))
                idx[i,j,k]=True

dep=-np.flip(m=ele)
res=np.array(object=[lon_res,lat_res,ele_res])
num=np.array(object=[lon_num,lat_num,ele_num])

vel[:,:,:]=np.flip(m=vel,axis=2) #(lon,lat,dep)
idx[:,:,:]=np.flip(m=idx,axis=2) #(lon,lat,dep)
np.savez(file=path_save,vel=vel,idx=idx,lon=lon,lat=lat,dep=dep,res=res,num=num)
