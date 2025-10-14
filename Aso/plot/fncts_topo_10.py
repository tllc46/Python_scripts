import numpy as np

topo=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/topo/Aso_topo_10m.npz")
topo=topo["topo"] #(6000,6750)=(lat,lon) [m]
grad=np.load(file="/home/tllc46/48NAS1/tllc46/Aso/topo/Aso_topo_10m_grad.npz")
grad=grad["grad"] #(6000,6750)=(lat,lon)
lat=np.arange(stop=6000)/9000+(32.5835-1/9000)
lon=np.arange(stop=6750)/9000+(130.7505-4/9000)

def init_prfl(coord_min,coord_max,coord_prfl):
    global lon,topo_prfl_lon
    global lat,topo_prfl_lat
    global topo,grad

    idx_lon_min=int(np.floor(9000*(coord_min[0]-130.7505)+4))
    idx_lon_max=int(np.ceil(9000*(coord_max[0]-130.7505)+4))+1
    idx_lon_prfl=round(9000*(coord_prfl[0]-130.7505)+4)
    lon=lon[idx_lon_min:idx_lon_max]

    idx_lat_min=int(np.floor(9000*(coord_min[1]-32.5835)+1))
    idx_lat_max=int(np.ceil(9000*(coord_max[1]-32.5835)+1))+1
    idx_lat_prfl=round(9000*(coord_prfl[1]-32.5835)+1)
    lat=lat[idx_lat_min:idx_lat_max]

    topo_prfl_lon=-0.001*topo[idx_lat_prfl,idx_lon_min:idx_lon_max]
    topo_prfl_lat=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_prfl]
    topo=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_min:idx_lon_max]
    grad=grad[idx_lat_min:idx_lat_max,idx_lon_min:idx_lon_max]

def init_lklhd(coord_min,coord_max,idx_loc,grid):
    global lon,topo_prfl_lon
    global lat,topo_prfl_lat
    global topo

    grid_lon=grid["lon"]
    grid_lat=grid["lat"]
    res_inv=grid["res_inv"] #(lon,lat,dep)

    idx_lon_min=int(np.floor(9000*(coord_min[0]-130.7505)+4))
    idx_lon_max=int(np.ceil(9000*(coord_max[0]-130.7505)+4))+1
    lon_inc=9000/res_inv[0]
    idx_lon_prfl=round(9000*(grid_lon[0]-130.7505)+4+lon_inc*idx_loc[0])
    lon=lon[idx_lon_min:idx_lon_max]

    idx_lat_min=int(np.floor(9000*(coord_min[1]-32.5835)+1))
    idx_lat_max=int(np.ceil(9000*(coord_max[1]-32.5835)+1))+1
    lat_inc=9000/res_inv[1]
    idx_lat_prfl=round(9000*(grid_lat[0]-32.5835)+1+lat_inc*idx_loc[1])
    lat=lat[idx_lat_min:idx_lat_max]

    topo_prfl_lon=-0.001*topo[idx_lat_prfl,idx_lon_min:idx_lon_max]
    topo_prfl_lat=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_prfl]
    topo=-0.001*topo[idx_lat_min:idx_lat_max,idx_lon_min:idx_lon_max]
