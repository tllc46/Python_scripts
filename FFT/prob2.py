import numpy as np
import matplotlib.pyplot as plt

A=1
tau=20
T=256
dt=1
N=int(T/dt)

omega_analy=np.arange(stop=int(6*np.pi/0.01)+1)*0.01-3*np.pi
omega_dft=np.arange(stop=3*N+1)*2*np.pi/N-3*np.pi

analy=2*A*np.sin(0.5*tau*omega_analy)/omega_analy

seq=np.arange(stop=N)
seq_rect=np.concatenate((seq[:int(0.5*tau)+1],seq[-int(0.5*tau):]))
complex_dft=A*sum(np.exp(-1j*2*np.pi*seq*seq_rect[:,np.newaxis]/N))
dft=complex_dft.real
dft=np.concatenate((dft[int(0.5*N):],dft[:int(0.5*N)])) #omega 범위가 -pi ~ pi가 되도록 재배열
dft=np.tile(A=dft,reps=3) #omega 범위가 -3*pi ~ 3*pi가 되도록 3회 반복
dft=np.append(arr=dft,values=dft[0]) #omega=-3*pi의 요소를 뒤에 추가

fig=plt.figure()
ax=fig.subplots()
ax.scatter(x=omega_dft,y=dft)
ax.plot(omega_analy,analy,color="r")
ax.set_xlim(left=-3*np.pi,right=3*np.pi)
ax.set_xlabel(xlabel="$\\omega$(rad/s)")
ax.set_title(label=f"$\\tau$={tau}, A={A}")
ax.legend(labels=["DFT","analytic"])
plt.show()
