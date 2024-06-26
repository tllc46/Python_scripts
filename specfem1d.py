from os.path import exists
from os import mkdir

import numpy as np
from scipy.special import gamma,roots_jacobi,eval_legendre

#상수
fixed_bc=False

nspec=500 #총 요소 개수
ngll=4 #GLL 점 개수
nglob=(ngll-1)*nspec+1 #총 GLL 점 개수(중복 미포함)
cfl=0.4 #Courant-Friedrichs-Lewy(CFL) stability value
nstep=1641 #시간 단계 개수
nsnap=200 #결과 출력 시간 단계

length=50E3 #[m]
density=1500 #[kg/m^3]
rigidity=30 #[Pa]
C=0 #C*v 감쇠

ispec_source=0 #source의 요소 번호
i_source=0 #source의 GLL 점 번호
hdur_divided_by_dt_of_the_source=82.035395797065604
source_amp=1E7
ireceiver=2*ngll-3 #지진계 위치 색인

#Newmark-beta method
#in case of explicit central difference scheme
newmark_beta=0
newmark_gamma=0.5

def source_time_function(t,hdur):
    alpha=2.628/hdur
    return -2*(alpha**3)*t*np.exp(-(alpha*t)**2)/np.sqrt(np.pi)

#중요: n은 GLL 점 개수(ngll)-1, 원래 코드와 통일
gll_points={2:np.array([-1,0,1]),
        3:np.array([-1,-np.sqrt(1/5),np.sqrt(1/5),1]),
        4:np.array([-1,-np.sqrt(3/7),0,np.sqrt(3/7),1]),
        5:np.array([-1,-np.sqrt(1/3*(1+np.sqrt(4/7))),-np.sqrt(1/3*(1-np.sqrt(4/7))),np.sqrt(1/3*(1-np.sqrt(4/7))),np.sqrt(1/3*(1+np.sqrt(4/7))),1]),
        6:np.array([-1,-np.sqrt(5/11*(1+np.sqrt(4/15))),-np.sqrt(5/11*(1-np.sqrt(4/15))),0,np.sqrt(5/11*(1-np.sqrt(4/15))),np.sqrt(5/11*(1+np.sqrt(4/15))),1]),
        7:np.array([-1,-0.8717401485,-0.5917001814,-0.2092992179,0.2092992179,0.5917001814,0.8717401485,1])}

#중요: n은 GLL 점 개수(ngll)-1, 원래 코드와 통일
gll_weights={2:np.array([1/3,4/3,1/3]),
        3:np.array([1/6,5/6,5/6,1/6]),
        4:np.array([1/10,49/90,32/45,49/90,1/10]),
        5:np.array([1/15,(14-np.sqrt(7))/30,(14+np.sqrt(7))/30,(14+np.sqrt(7))/30,(14-np.sqrt(7))/30,1/15]),
        6:np.array([1/21,(124-7*np.sqrt(15))/350,(124+7*np.sqrt(15))/350,256/525,(124+7*np.sqrt(15))/350,(124-7*np.sqrt(15))/350,1/21]),
        7:np.array([1/28,0.2107042271,0.3411226924,0.4124587946,0.4124587946,0.3411226924,0.2107042271,1/28])}

#중요: n은 GLL 점 개수(ngll)-1, 원래 코드와 통일
lagrange_deriv={2:np.array(
        [[-3/2,2,-1/2],
         [-1/2,0,1/2],
         [1/2,-2,3/2]]),
        3:np.array(
        [[-3,5*(np.sqrt(5)+1)/4,-5*(np.sqrt(5)-1)/4,1/2],
         [-(np.sqrt(5)+1)/4,0,np.sqrt(5)/2,-(np.sqrt(5)-1)/4],
         [(np.sqrt(5)-1)/4,-np.sqrt(5)/2,0,(np.sqrt(5)+1)/4],
         [-1/2,5*(np.sqrt(5)-1)/4,-5*(np.sqrt(5)+1)/4,3]])}

def endw1(n,alpha,beta): #xi=-1 GLL 점의 가중치
    #중요: n은 GLL 점 개수(ngll)-1, 원래 코드와 통일
    return 2**(alpha+beta+1)*(beta+1)*gamma(z=beta+1)**2*gamma(z=n)*gamma(z=n+alpha+1)/(gamma(z=n+beta+1)*gamma(z=n+alpha+beta+2))

def endw2(n,alpha,beta): #xi=1 GLL 점의 가중치
    #중요: n은 GLL 점 개수(ngll)-1, 원래 코드와 통일
    return 2**(alpha+beta+1)*(alpha+1)*gamma(z=alpha+1)**2*gamma(z=n)*gamma(z=n+beta+1)/(gamma(z=n+alpha+1)*gamma(z=n+alpha+beta+2))

#zwgljd
xigll,wgll=roots_jacobi(n=ngll-2,alpha=1,beta=1) #xi_gll, weight_gll
wgll/=1-xigll**2
xigll=np.insert(arr=xigll,obj=[0,ngll-2],values=[-1,1]) #GLL 점에 -1, 1도 추가
wgll=np.insert(arr=wgll,obj=[0,ngll-2],values=[endw1(n=ngll-1,alpha=0,beta=0),endw2(n=ngll-1,alpha=0,beta=0)]) #가중치에 -1,1의 가중치도 추가

#defineDerivationMatrices
hprime=[]
index=0
while index<int(ngll**2/2): #hprime은 부호만 반대일 뿐 점대칭 행렬이므로 반만 계산
    i=int(index/ngll)
    j=index%ngll
    if i==j==0:
        hprime.append(-0.25*ngll*(ngll-1))
    elif i==j and 0<i and 0<j:
        hprime.append(0)
    else:
        hprime.append(eval_legendre(n=ngll-1,x=xigll[i])/(eval_legendre(n=ngll-1,x=xigll[j])*(xigll[i]-xigll[j]))) #l_j'(xi_i)
    index+=1

if ngll%2: #ngll이 홀수
    hprime.append(0)

hprime=np.array(object=hprime)
hprime=np.append(arr=hprime,values=-hprime[int(ngll**2/2)-1::-1]).reshape(ngll,ngll) #중요: hprime(i,j)=l_j'(xi_i), 차원은 (i,j)

#makeGrid
anchor=np.arange(stop=nspec+1)*length/nspec #요소 경계 좌표
rho=np.ones(shape=(ngll,nspec))*density #밀도
mu=np.ones(shape=(ngll,nspec))*rigidity #강성 계수
dxi_dx=2/(anchor[1:]-anchor[:-1]) #dxi/dx
x=(0.5*(anchor[1:]+anchor[:-1]))[:,np.newaxis] + (0.5*(anchor[1:]-anchor[:-1]))[:,np.newaxis]*xigll[:-1]
x=np.append(arr=x.flatten(),values=length) #GLL 점 좌표

#estimateTimeStep
dh=length/(nglob-1)
v=np.sqrt(rigidity/density)
deltat=cfl*dh/v
hdur=hdur_divided_by_dt_of_the_source*deltat

#assembleGlobalMassMatrix
mass_global=np.zeros(shape=nglob)
for i in range(nspec):
    mass_local=wgll*(rho[:,i]+0.5*C*deltat)/dxi_dx[i]
    mass_global[(ngll-1)*i:(ngll-1)*i+ngll]+=mass_local

#initializeTimeMarching
displ=np.zeros(shape=nglob)
veloc=np.zeros(shape=nglob)
accel=np.zeros(shape=nglob)
if not exists(path=f"OUTPUT_FILES"):
    mkdir(path=f"OUTPUT_FILES")

seismogram=[]
#mainTimeLoop
for it in range(nstep):
    displ+=deltat*veloc+0.5*deltat**2*(1-2*newmark_beta)*accel
    veloc+=(1-newmark_gamma)*deltat*accel
    accel=np.zeros(shape=nglob)
    for i in range(nspec):
        temp=np.matmul(displ[(ngll-1)*i:(ngll-1)*i+ngll],hprime.T)*wgll*mu[:,i]*dxi_dx[i] #temp(dim=(1,k)) = u(xi_i)(dim=(1,i)) * l_i'(xi_k)(dim=(i,k)), Comput. Seis. p. 195
        templ=np.matmul(temp,hprime) #templ(dim=(1,j)) = temp(dim=(1,k)) * l_j'(xi_k)(dim=(k,j)), Comput. Seis. p. 195
        accel[(ngll-1)*i:(ngll-1)*i+ngll]-=templ

    #addSourceAtGlobalLevel
    accel[(ngll-1)*ispec_source+i_source]+=source_amp*source_time_function(t=it*deltat-hdur,hdur=hdur)
    if fixed_bc:
        accel[0]=0
        accel[-1]=0
    accel-=C*veloc
    accel/=mass_global
    veloc+=newmark_gamma*deltat*accel

    #writeOutSnapshots
    if not (it+1)%nsnap:
        with open(file=f"OUTPUT_FILES/snapshot_forward_normal{it+1:05}",mode="w") as file:
            for i in range(nglob):
                file.write(f"{x[i]} {displ[i]}\n")

    seismogram.append(displ[ireceiver])
