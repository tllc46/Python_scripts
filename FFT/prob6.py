import numpy as np
import matplotlib.pyplot as plt

a_n=np.array(object=[2,1.5,1])
omega_n=np.array(object=[np.pi/6,np.pi/3,0.5*np.pi])
phi_n=np.array(object=[0,np.pi/6,np.pi/3])
t_n=np.array(object=[90,130,170])
T=256
dt=0.25
dt_noise=0.5
N=int(T/dt)
N_noise=int(T/dt_noise)

#6-1
seq=np.arange(stop=N)
time=dt*seq
x=sum(a_n[:,np.newaxis]*np.cos(omega_n[:,np.newaxis]*time+phi_n[:,np.newaxis])*np.exp(-((time-t_n[:,np.newaxis])/35)**2))
noise=np.zeros(shape=N)
noise[::2]=np.random.uniform(low=-0.2,high=0.2,size=N_noise)
x_noise=x+noise

omega=np.arange(stop=N+1)*8*np.pi/N-4*np.pi

#6-2
complex_dft=sum(x_noise[:,np.newaxis]*np.exp(-1j*2*np.pi*seq*seq[:,np.newaxis]/N))
abs_dft=0.25*abs(complex_dft) #normalization
abs_dft=np.concatenate((abs_dft[int(0.5*N):],abs_dft[:int(0.5*N)+1])) #rearrange to make omega's range -4*pi~4*pi

#6-3
lp=np.zeros(shape=N)
ind_lp=int(np.floor((0.5*np.pi+3*np.sqrt(2)/35)/(2*np.pi/256))) #저역 통과 기준 진동수: 0.5*pi+3*(표준 편차)
lp[:ind_lp+1]=1
lp[-ind_lp:]=1
lp_dft=complex_dft*lp
x_lp=sum(lp_dft[:,np.newaxis]*np.exp(1j*2*np.pi*seq*seq[:,np.newaxis]/N)).real/N

#6-4
hp=np.zeros(shape=N)
hp[ind_lp+1:-ind_lp]=1
hp_dft=complex_dft*hp
x_hp=sum(hp_dft[:,np.newaxis]*np.exp(1j*2*np.pi*seq*seq[:,np.newaxis]/N)).real/N

#6-5
bp=np.zeros(shape=N)
ind_low=int(np.ceil(0.25*np.pi/(2*np.pi/256)))
ind_high=int(np.floor((2*np.pi/3)/(2*np.pi/256)))
bp[ind_low:ind_high+1]=1
bp[-ind_high:-ind_low+1]=1
bp_dft=complex_dft*bp
x_bp=sum(bp_dft[:,np.newaxis]*np.exp(1j*2*np.pi*seq*seq[:,np.newaxis]/N)).real/N

#6-6
omega_sym=np.concatenate((omega[int(0.5*N):-1],omega[int(0.5*N)+1:][::-1]))
gauss=np.exp(-0.25*(35*(omega_sym-omega_n[:,np.newaxis]))**2)
gauss_dft=complex_dft*gauss
x_gauss=np.sum(a=gauss_dft[:,:,np.newaxis]*np.exp(1j*2*np.pi*seq*seq[:,np.newaxis]/N),axis=1).real/N

fig=plt.figure()
ax=fig.subplots()
ax.plot(time,x_noise)
ax.set_xlim(left=0,right=256)
ax.set_xlabel(xlabel="time(s)")
ax.set_title(label="original")
plt.show()

fig=plt.figure()
ax=fig.subplots()
ax.plot(omega,abs_dft)
ax.set_xlim(left=-4*np.pi,right=4*np.pi)
ax.set_xlabel(xlabel="$\\omega$(rad/s)")
ax.set_title(label="amplitude spectrum")
plt.show()

fig=plt.figure()
gs=fig.add_gridspec(nrows=4,ncols=1)
ax1=fig.add_subplot(gs[:3,0])
ax1.plot(omega,abs_dft)
ax2=fig.add_subplot(gs[3,0],sharex=ax1)
ax2.plot(omega,np.concatenate((lp[int(0.5*N):],lp[:int(0.5*N)+1])),color="r")
ax2.set_xlim(left=-4*np.pi,right=4*np.pi)
ax2.set_xlabel(xlabel="$\\omega$(rad/s)")
ax1.label_outer()
plt.show()

fig=plt.figure()
ax=fig.subplots()
ax.plot(time,x_lp)
ax.plot(time,x)
ax.set_xlim(left=0,right=256)
ax.set_xlabel(xlabel="time(s)")
ax.set_title(label="low pass filter")
ax.legend(labels=["filtered","theoretical"])
plt.show()

fig=plt.figure()
gs=fig.add_gridspec(nrows=4,ncols=1)
ax1=fig.add_subplot(gs[:3,0])
ax1.plot(omega,abs_dft)
ax2=fig.add_subplot(gs[3,0],sharex=ax1)
ax2.plot(omega,np.concatenate((hp[int(0.5*N):],hp[:int(0.5*N)+1])),color="r")
ax2.set_xlim(left=-4*np.pi,right=4*np.pi)
ax2.set_xlabel(xlabel="$\\omega$(rad/s)")
ax1.label_outer()
plt.show()

fig=plt.figure()
ax=fig.subplots()
ax.plot(time,x_hp)
ax.set_xlim(left=0,right=256)
ax.set_xlabel(xlabel="time(s)")
ax.set_title(label="high pass filter")
plt.show()

fig=plt.figure()
gs=fig.add_gridspec(nrows=4,ncols=1)
ax1=fig.add_subplot(gs[:3,0])
ax1.plot(omega,abs_dft)
ax2=fig.add_subplot(gs[3,0],sharex=ax1)
ax2.plot(omega,np.concatenate((bp[int(0.5*N):],bp[:int(0.5*N)+1])),color="r")
ax2.set_xlim(left=-4*np.pi,right=4*np.pi)
ax2.set_xlabel(xlabel="$\\omega$(rad/s)")
ax1.label_outer()
plt.show()

fig=plt.figure()
ax=fig.subplots()
ax.plot(time,x_bp)
ax.set_xlim(left=0,right=256)
ax.set_xlabel(xlabel="time(s)")
ax.set_title(label="band pass filter")
plt.show()

fig=plt.figure()
gs=fig.add_gridspec(nrows=6,ncols=1)
ax1=fig.add_subplot(gs[:3,0])
ax1.plot(omega,abs_dft)
ax1.label_outer()
for i in range(3,6):
    ax=fig.add_subplot(gs[i,0],sharex=ax1)
    ax.plot(omega,np.concatenate((gauss[i-3,int(0.5*N):],gauss[i-3,:int(0.5*N)+1])),color="r")
    if i!=5:
        ax.label_outer()
    else:
        ax.set_xlim(left=-4*np.pi,right=4*np.pi)
        ax.set_xlabel(xlabel="$\\omega$(rad/s)")
plt.show()

fig=plt.figure()
ax_array=fig.subplots(nrows=3,ncols=1,sharex=True)
for i in range(3):
    ax_array[i].plot(time,a_n[i]*np.cos(omega_n[i]*time+phi_n[i])*np.exp(-((time-t_n[i])/35)**2))
    ax_array[i].plot(time,x_gauss[i])
ax_array[0].set_title(label="gaussian filter")
ax_array[0].legend(labels=["theoretical","filtered"])
ax_array[2].set_xlim(left=0,right=256)
ax_array[2].set_xlabel(xlabel="time(s)")
plt.show()
