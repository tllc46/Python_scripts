import numpy as np
import matplotlib.pyplot as plt

A=1
tau=20
T=256
dt=1
N=int(T/dt)

omega_analy=np.arange(start=-3*np.pi,stop=3*np.pi,step=0.01)
omega_dft=3*np.pi*np.linspace(start=-1,stop=1,num=3*N+1)

analy=2*A*np.sin(0.5*tau*omega_analy)/omega_analy

seq=np.arange(stop=N)
seq_rect=np.concatenate((seq[:int(0.5*tau)+1],seq[-int(0.5*tau):]))
complex_dft=A*sum(np.exp(-1j*2*np.pi*seq*seq_rect[:,np.newaxis]/N))
dft=complex_dft.real
dft=np.concatenate((dft[int(0.5*N):],dft[:int(0.5*N)])) #rearange to make omega's range -pi~pi
dft=np.tile(A=dft,reps=3) #repeat 3 times to make omega's range -3*pi~3*pi
dft=np.append(arr=dft,values=dft[0]) #append first element to complete end of omega's range

fig=plt.figure()
ax=fig.subplots()
ax.scatter(x=omega_dft,y=dft)
ax.plot(omega_analy,analy,color="r")
ax.set_xlim(left=-3*np.pi,right=3*np.pi)
ax.set_xlabel(xlabel="$\omega$(rad/s)")
ax.set_title(label=f"$\\tau$={tau}, A={A}")
ax.legend(labels=["DFT","analytic"])
plt.show()
