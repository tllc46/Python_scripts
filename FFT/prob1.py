import numpy as np
import matplotlib.pyplot as plt

A=1
omega0=0.5
T=128
dt=1
N=T//dt

def analytic(omega):
    if omega==omega0 or omega==-omega0:
        return A*np.pi
    else:
        return 0

omega_analy=np.arange(stop=int(6/omega0)+1)*omega0-3
omega_dft=np.arange(stop=N+1)*2*np.pi/N-np.pi

abs_analy=[analytic(omega=omega) for omega in omega_analy]
angle_analy=np.zeros(shape=len(omega_analy))

seq=np.arange(stop=N)
complex_dft=sum(A*np.cos(omega0*seq[:,np.newaxis])*np.exp(-1j*2*np.pi*seq*seq[:,np.newaxis]/N))
abs_dft=2*np.pi/128*abs(complex_dft) #정규화
angle_dft=np.angle(z=complex_dft)
abs_dft=np.concatenate((abs_dft[N//2:],abs_dft[:N//2+1])) #omega 범위가 -pi ~ pi가 되도록 재배열
angle_dft=np.concatenate((angle_dft[N//2:],angle_dft[:N//2+1]))

fig=plt.figure(1)
ax=fig.subplots()
ax.axvline(x=-0.5,linestyle=":",zorder=1)
ax.axvline(x=0.5,linestyle=":",zorder=1)
pathcol_dft=ax.scatter(x=omega_dft,y=abs_dft)
pathcol_analy=ax.scatter(x=omega_analy,y=abs_analy)
ax.set_xlim(left=-np.pi,right=np.pi)
ax.set_xlabel(xlabel="$\\omega$(rad/s)")
ax.set_title(label="absolute(magnitude, modulus)")
ax.legend(handles=[pathcol_dft,pathcol_analy],labels=["DFT","analytic"])
plt.show()

fig=plt.figure(2)
ax=fig.subplots()
ax.scatter(x=omega_dft,y=angle_dft)
ax.axvline(x=-0.5,linestyle=":",zorder=1)
ax.axvline(x=0.5,linestyle=":",zorder=1)
ax.scatter(x=omega_analy,y=angle_analy)
ax.set_xlim(left=-np.pi,right=np.pi)
ax.set_xlabel(xlabel="$\\omega$(rad/s)")
ax.set_title(label="angle(argument)(rad)")
ax.legend(labels=["DFT","analytic"])
plt.show()
