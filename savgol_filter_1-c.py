#implementation of scipy.signal.savgol_filter() when polyorder=1
#equivalent to y=scipy.signal.savgol_filter(x=data,window_length=window_length,polyorder=polyorder)
#savgol filter reduces to moving average

import numpy as np

window_length=11
polyorder=1
data=np.array(object=[1,4,3,5,8,3,2,4,5,8,9,3,5,3,2,1,4,2,6,8,3])

npts=len(data)
halflen=window_length//2

y=np.empty(shape=npts)

#moving average
data_cumsum=np.cumsum(a=data)
y[halflen+1:-halflen]=data_cumsum[window_length:]-data_cumsum[:-window_length]
y[halflen]=data_cumsum[window_length-1]
y[halflen:-halflen]/=window_length

#_fit_edge
#in case of polyorder=1, just linear regression
t_half=np.arange(start=-halflen,stop=halflen+1)
denom=halflen*(halflen+1)*window_length//3

#left
slope=sum(t_half*data[:window_length])/denom
y[:halflen]=np.mean(a=data[:window_length])+slope*t_half[:halflen]

#right
slope=sum(t_half*data[-window_length:])/denom
y[-halflen:]=np.mean(a=data[-window_length:])+slope*t_half[-halflen:]
