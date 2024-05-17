import numpy as np
import matplotlib.pyplot as plt

A=1
omega0=0.5
T=128
dt=1
N=int(T/dt)

def analytic(omega):
    if omega==omega0 or omega==-omega0:
        return A*np.pi
    else:
        return 0

omega_analy=np.arange(start=-3,stop=3+omega0,step=omega0)
omega_dft=np.pi*np.linspace(start=-1,stop=1,num=N+1)

abs_analy=[analytic(omega=omega) for omega in omega_analy]
angle_analy=np.zeros(shape=len(omega_analy))

seq=np.arange(stop=N)
complex_dft=sum(A*np.cos(omega0*seq[:,np.newaxis])*np.exp(-1j*2*np.pi*seq*seq[:,np.newaxis]/N))
abs_dft=2*np.pi/128*abs(complex_dft) #normalization
angle_dft=np.angle(z=complex_dft)
abs_dft=np.concatenate((abs_dft[int(0.5*N):],abs_dft[:int(0.5*N)+1])) #rearrange to make omega's range -pi~pi
angle_dft=np.concatenate((angle_dft[int(0.5*N):],angle_dft[:int(0.5*N)+1]))

fig=plt.figure(1)
ax=fig.subplots()
ax.axvline(x=-0.5,linestyle=":",zorder=1)
ax.axvline(x=0.5,linestyle=":",zorder=1)
pathcol_dft=ax.scatter(x=omega_dft,y=abs_dft)
pathcol_analy=ax.scatter(x=omega_analy,y=abs_analy)
ax.set_xlim(left=-np.pi,right=np.pi)
ax.set_xlabel(xlabel="$\omega$(rad/s)")
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
ax.set_xlabel(xlabel="$\omega$(rad/s)")
ax.set_title(label="angle(argument)(rad)")
ax.legend(labels=["DFT","analytic"])
plt.show()
