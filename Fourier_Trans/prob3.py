import numpy as np
import matplotlib.pyplot as plt

A=10
tau=2
dT=20
dt=1
N=int(dT/dt)

def analytic(omega):
    if omega:
        return 2*A*np.sin(0.5*tau*omega)/omega
    else:
        return A*tau

omega=np.arange(stop=N+1)*2*np.pi/N-np.pi

dirac=20*2*np.pi/dT*np.ones(shape=len(omega))

analy=np.array(object=[analytic(omega=omega_) for omega_ in omega])
analy*=0.1*np.pi #정규화

seq=np.arange(stop=N)
seq_rect=np.concatenate((seq[:int(0.5*tau)+1],seq[-int(0.5*tau):]))
complex_dft=A*sum(np.exp(-1j*2*np.pi*seq*seq_rect[:,np.newaxis]/N))
dft=0.1*np.pi*complex_dft.real #normalization
dft=np.concatenate((dft[int(0.5*N):],dft[:int(0.5*N)+1])) #omega 범위가 -pi ~ pi가 되도록 재배열

fig=plt.figure()
ax=fig.subplots()
ax.scatter(x=omega,y=dirac)
ax.scatter(x=omega,y=analy)
ax.scatter(x=omega,y=dft)
ax.set_xlim(left=-np.pi,right=np.pi)
ax.set_xlabel(xlabel="$\omega$(rad/s)")
ax.legend(labels=["Dirac $\delta$","analytic","DFT"])
plt.show()
