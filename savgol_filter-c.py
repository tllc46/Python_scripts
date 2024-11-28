#implementation of scipy.signal.savgol_filter()

import numpy as np
from numpy.linalg import inv
from scipy.linalg import lstsq

window_length=11
polyorder=1
data=np.array(object=[1,4,3,5,8,3,2,4,5,8,9,3,5,3,2,1,4,2,6,8,3])

npts=len(data)
halflen=window_length//2

def savgol_coeffs():
    rem=window_length%2

    if rem:
        pos=halflen
    else:
        pos=halflen-0.5

    x=pos-np.arange(stop=window_length)
    order=np.arange(stop=polyorder+1)[:,None]
    A=x**order #(polyorder+1,window_length)

    y=np.zeros(shape=polyorder+1)
    y[0]=1
    coeffs,_,_,_=lstsq(a=A,b=y)

    return coeffs

coeffs=savgol_coeffs()

#convolve_1d
y=np.empty(shape=npts)
for i in range(halflen,npts-halflen):
    start_idx=i-(window_length-1)//2
    y[i]=sum(data[start_idx:start_idx+window_length]*coeffs)

#_fit_edge
x=np.arange(stop=window_length)[:,None]
order=np.arange(stop=polyorder+1)
A=x**order #(window_length,polyorder+1)

sol_mat=np.matmul(inv(a=np.matmul(A.T,A)),A.T) #(polyorder+1,window_length)
left=x[:halflen]
right=x[-halflen:]

left=left**order #(halflen,polyorder+1)
right=right**order #(halflen,polyorder+1)

left=np.matmul(left,sol_mat) #(halflen,window_length)
right=np.matmul(right,sol_mat) #(halflen,window_length)

left_edge=np.matmul(left,data[:window_length]) #(halflen,)
right_edge=np.matmul(right,data[-window_length:]) #(halflen,)

y[:halflen]=left_edge
y[-halflen:]=right_edge
