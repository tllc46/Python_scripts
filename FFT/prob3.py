import numpy as np
import matplotlib.pyplot as plt

A=10
tau=2
dT=20
dt=1
N=dT//dt

def analytic(omega):
    if omega:
        return 2*A*np.sin(0.5*tau*omega)/omega
    else:
        return A*tau

omega=np.arange(stop=N+1)*2*np.pi/N-np.pi

dirac=A*tau*2*np.pi/dT*np.ones(shape=len(omega))

analy=np.array(object=[analytic(omega=omega_) for omega_ in omega])
analy*=0.1*np.pi #정규화

seq=np.arange(stop=N)
if int(tau/2)<1:
    seq_rect=np.array([0])
else:
    seq_rect=np.concatenate((seq[:int(tau/2)+1],seq[-int(tau/2):]))complex_dft=A*sum(np.exp(-1j*2*np.pi*seq*seq_rect[:,np.newaxis]/N))
dft=0.1*np.pi*complex_dft.real #정규화
dft=np.concatenate((dft[N//2:],dft[:N//2+1])) #omega 범위가 -pi ~ pi가 되도록 재배열

fig=plt.figure()
ax=fig.subplots()
ax.scatter(x=omega,y=dirac)
ax.scatter(x=omega,y=analy)
ax.scatter(x=omega,y=dft)
ax.set_xlim(left=-np.pi,right=np.pi)
ax.set_xlabel(xlabel="$\\omega$(rad/s)")
ax.legend(labels=["Dirac $\\delta$","analytic","DFT"])
plt.show()
