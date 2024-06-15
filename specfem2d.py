from os.path import exists
from os import mkdir

import numpy as np
from scipy.special import gamma,roots_jacobi,eval_legendre

#constants
rho=2700 #density [kg/m^3]
cp=3000 #P wave velocity [m/s]
cs=cp/np.sqrt(3) #S wave velocity [m/s]

xmin=0 #left [m]
xmax=4000 #right [m]
zmin=0 #bottom [m]
zmax=3000 #top [m]
nelem_x=80 #no. of spectral elements along x
nelem_z=60 #no. of spectral elements along z
nspec=nelem_x*nelem_z #no. of total spectral elements
npgeo=(nelem_x+1)*(nelem_z+1) #no. of total nodes

ngllx=5 #no. of gll points along x
ngllz=5 #no. of gll points along z
nxyz=ngllx*ngllz #total no. of gll points in 1 spectral element
ntot=nxyz*nspec #total no. of gll points (includes overlap)
nglob=(nelem_x*(ngllx-1)+1)*(nelem_z*(ngllz-1)+1) #no. of unique gll points

ngnod=4 #need 4 nodes to describe spectral element mesh

deltat=0.0011 #[s]
nstep=1600 #no. of time steps
nstep_between_output_info=100

#Newmark-beta method
#in case of explicit central difference scheme
newmark_beta=0
newmark_gamma=0.5

nspec_source=int(nspec/2)-int(nelem_x/2)-1
igll_source=ngllx-1
jgll_source=ngllz-1

nspec_receiver=int(2*nspec/3)-int(nelem_x/4)-1
igll_receiver=0
jgll_receiver=0

f0=10
t0=1.2/f0
factor_amplitude=1E+10
a=(np.pi*f0)**2

ndim=2
smallvaltol=10
stability_threshold=1E+25

knods=np.empty(shape=(ngnod,nspec),dtype=int)
coorg=np.empty(shape=(ndim,npgeo))
shape2D=np.empty(shape=(ngnod,ngllx,ngllz))
dershape2D_z=np.empty(shape=(ngnod,ngllz))
dershape2D_x=np.empty(shape=(ngnod,ngllx))
xp=np.empty(shape=(ngllx,ngllz,nspec))
zp=np.empty(shape=(ngllx,ngllz,nspec))
ibool=np.empty(shape=ntot,dtype=int)
rmass_inverse=np.zeros(shape=nglob)
displ=np.zeros(shape=(ndim,nglob))
veloc=np.zeros(shape=(ndim,nglob))
accel=np.zeros(shape=(ndim,nglob))
seismogram=[]

def endw1(n,alpha,beta): #xi=-1 end point weight
    #important: n must be no. of gll points (ngll) - 1, consistent with Fortran code
    return 2**(alpha+beta+1)*(beta+1)*gamma(z=beta+1)**2*gamma(z=n)*gamma(z=n+alpha+1)/(gamma(z=n+beta+1)*gamma(z=n+alpha+beta+2))

def endw2(n,alpha,beta): #xi=1 end point weight
    #important: n must be no. of gll points (ngll) - 1, consistent with Fortran code
    return 2**(alpha+beta+1)*(alpha+1)*gamma(z=alpha+1)**2*gamma(z=n)*gamma(z=n+beta+1)/(gamma(z=n+alpha+1)*gamma(z=n+alpha+beta+2))

def zwgljd(n,alpha,beta): #zeros and weights of glj points
    z,w=roots_jacobi(n=n-2,alpha=alpha+1,beta=beta+1) #zeros, weights
    w/=1-z**2
    z=np.insert(arr=z,obj=[0,n-2],values=[-1,1]) #add -1 and 1
    w=np.insert(arr=w,obj=[0,n-2],values=[endw1(n=n-1,alpha=alpha,beta=beta),endw2(n=n-1,alpha=alpha,beta=beta)]) #add -1 and 1's weights
    return z,w

def define_derivation_matrices(zgll,ngll):
    #zgll: zeros (gll points)
    #ngll: no. of gll points
    hprime=[] #h'
    index=0
    while index<int(ngll**2/2): #because hprime is symmetric with just opposite sign, calculate only half
        i=int(index/ngll) #row no.
        j=index%ngll #column no.
        #lagrange_deriv_GLL
        if i==j==0:
            hprime.append(-0.25*ngll*(ngll-1))
        elif i==j and 0<i and 0<j:
            hprime.append(0)
        else:
            hprime.append(eval_legendre(n=ngll-1,x=zgll[i])/(eval_legendre(n=ngll-1,x=zgll[j])*(zgll[i]-zgll[j]))) #l_j'(xi_i)
        index+=1

    if ngll%2: #ngll is odd
        hprime.append(0)

    hprime=np.array(object=hprime)
    hprime=np.append(arr=hprime,values=-hprime[int(ngll**2/2)-1::-1]).reshape(ngll,ngll) #important: hprime(i,j)=l_j'(xi_i), dim=i x j
    return hprime

xigll,wxgll=zwgljd(n=ngllx,alpha=0,beta=0)
zigll,wzgll=zwgljd(n=ngllz,alpha=0,beta=0)
hprime_xx=define_derivation_matrices(zgll=xigll,ngll=ngllx)
hprime_zz=define_derivation_matrices(zgll=zigll,ngll=ngllz)
hprimewgll_xx=hprime_xx*wxgll[:,np.newaxis] #(ngllx,ngllx) * (ngllx,1) = (ngllx,ngllx)
hprimewgll_zz=hprime_zz*wzgll[:,np.newaxis] #(ngllz,ngllz) * (ngllz,1) = (ngllz,ngllz)

#important: global node numbering,
#starts at (xmin,zmin), increases x direction, consistent with Fortran code
nod_num=np.arange(stop=npgeo).reshape(nelem_z+1,nelem_x+1)

#each spectral element's global node numbers
#important: spectral element's order,
#starts at bottom left, increases x direction, consistent with Fortran code
#important: in case of 4 nodes,
#xi coordinate (-1,-1) -> (1,-1) -> (1,1) -> (-1,1) order, consistent with Fortran code
knods[0]=nod_num[:-1,:-1].flatten() #(-1,-1)
knods[1]=nod_num[:-1,1:].flatten() #(1,-1)
knods[2]=nod_num[1:,1:].flatten() #(1,1)
knods[3]=nod_num[1:,:-1].flatten() #(-1,1)

#mapping global node number to node's real coordinate
xrow=xmin+np.arange(stop=nelem_x+1)*(xmax-xmin)/nelem_x
coorg[0]=np.tile(A=xrow,reps=nelem_z+1) #x coordinate
zcol=zmin+np.arange(stop=nelem_z+1)*(zmax-zmin)/nelem_z
coorg[1]=np.repeat(a=zcol,repeats=nelem_x+1) #z coordinate

#define_shape_functions
shape2D[0]=0.25*(xigll[:,np.newaxis]-1)*(zigll-1) #(-1,-1), (ngllx,1) + (ngllz) = (ngllx,ngllz)
shape2D[1]=-0.25*(xigll[:,np.newaxis]+1)*(zigll-1) #(1,-1)
shape2D[2]=0.25*(xigll[:,np.newaxis]+1)*(zigll+1) #(1,1)
shape2D[3]=-0.25*(xigll[:,np.newaxis]-1)*(zigll+1) #(-1,1)

if (sum(shape2D)!=np.ones(shape=(5,5))).any():
    print("sum of shape functions are not 1")
    exit()

dershape2D_z[0]=0.25*(zigll-1)
dershape2D_z[1]=-0.25*(zigll-1)
dershape2D_z[2]=0.25*(zigll+1)
dershape2D_z[3]=-0.25*(zigll+1)

if any(sum(dershape2D_z)):
    print("sum of derivative gamma shape functions are not 0")
    exit()

dershape2D_x[0]=0.25*(xigll-1)
dershape2D_x[1]=-0.25*(xigll+1)
dershape2D_x[2]=0.25*(xigll+1)
dershape2D_x[3]=-0.25*(xigll-1)

if any(sum(dershape2D_x)):
    print("sum of derivative xi shape functions are not 0")

temp_xp=shape2D[:,:,:,np.newaxis]*coorg[0,knods][:,np.newaxis,np.newaxis] #(ngnod,ngllx,ngllz,1) * (ngnod,1,1,nspec) = (ngnod,ngllx,ngllz,nspec)
xp=sum(temp_xp) #(ngllx,ngllz,nspec)
temp_zp=shape2D[:,:,:,np.newaxis]*coorg[1,knods][:,np.newaxis,np.newaxis]
zp=sum(temp_zp)

x_source=xp[igll_source,jgll_source,nspec_source]
z_source=zp[igll_source,jgll_source,nspec_source]

x_receiver=xp[igll_receiver,jgll_receiver,nspec_receiver]
z_receiver=zp[igll_receiver,jgll_receiver,nspec_receiver]

xp=np.swapaxes(a=xp,axis1=0,axis2=2).flatten() #(nspec,ngllz,ngllx).flatten()
zp=np.swapaxes(a=zp,axis1=0,axis2=2).flatten()

xtypdist=min(max(xp)-min(xp),max(zp)-min(zp))
decimals=smallvaltol-int(np.log10(xtypdist))
xp=np.round(a=xp,decimals=decimals)
zp=np.round(a=zp,decimals=decimals)

#sort lexicographically, first by x coordinates, then by z to break tie
locval=np.lexsort(keys=(zp,xp))
xp=xp[locval]
zp=zp[locval]

ig=0
ibool[locval[0]]=ig
for i in range(1,ntot):
    if xp[i-1]<xp[i] or zp[i-1]<zp[i]:
        ig+=1
    ibool[locval[i]]=ig

if ig!=nglob-1:
    print("no. of unique gll points inconsistent")
    exit()

ibool=np.reshape(a=ibool,newshape=(nspec,ngllz,ngllx))
ibool=np.swapaxes(a=ibool,axis1=0,axis2=2) #(ngllx,ngllz,nspec)

#recompute_jacobian
temp_xxi=dershape2D_z[:,:,np.newaxis]*coorg[0,knods][:,np.newaxis] #(ngnod,ngllz,1) * (ngnod,1,nspec) = (ngnod,ngllz,nspec)
xxi=sum(temp_xxi) #x/xi, (ngllz,nspec)

temp_xxi=dershape2D_z[:,:,np.newaxis]*coorg[1,knods][:,np.newaxis]
zxi=sum(temp_xxi) #z/xi

temp_xgamma=dershape2D_x[:,:,np.newaxis]*coorg[0,knods][:,np.newaxis] #(ngnod,ngllx,1) * (ngnod,1,nspec) = (ngnod,ngllx,nspec)
xgamma=sum(temp_xgamma) #x/gamma, (ngllx,nspec)

temp_xgamma=dershape2D_x[:,:,np.newaxis]*coorg[1,knods][:,np.newaxis]
zgamma=sum(temp_xgamma) #z/gamma

jacobian=xxi*zgamma[:,np.newaxis]-xgamma[:,np.newaxis]*zxi #(ngllz,nspec) * (ngllx,1,nspec) = (ngllx,ngllz,nspec)

if (jacobian<=0).any():
    print("jacobian is less than or equal to 0")
    exit()

xix=zgamma[:,np.newaxis]/jacobian #xi/x, (ngllx,1,nspec) / (ngllx,ngllz,nspec) = (ngllx,ngllz,nspec)
gammax=-zxi/jacobian #gamma/x, (ngllz,nspec) / (ngllx,ngllz,nspec) = (ngllx,ngllz,nspec)
xiz=-xgamma[:,np.newaxis]/jacobian #xi/z
gammaz=xxi/jacobian #gamma/z

for i in range(nspec):
    rmass_inverse[ibool[:,:,i]]+=wxgll[:,np.newaxis]*wzgll*rho*jacobian[:,:,i] #(ngllx,1) * (ngllz) * (ngllx,ngllz) = (ngllx,ngllz)

rmass_inverse=1/rmass_inverse

mu=rho*cs**2
lambda_=rho*cp**2-2*mu
lambdaplus2mu=lambda_+2*mu

for it in range(nstep):
    if not (it+1)%nstep_between_output_info or it==4 or it==nstep-1:
        usolidnorm=max(np.sqrt(displ[0]**2+displ[1]**2))
        if usolidnorm<0 or stability_threshold<usolidnorm:
            print("maximum norm displacement is larger than stability threshold, code becomes unstable and blow up")
            exit()

    displ+=deltat*veloc+0.5*deltat**2*(1-2*newmark_beta)*accel
    veloc+=(1-newmark_gamma)*deltat*accel
    accel=np.zeros(shape=(ndim,nglob))
    for i in range(nspec):
        temp_dux_dxi=displ[0,ibool[:,:,i]][:,np.newaxis]*np.swapaxes(a=hprime_xx,axis1=0,axis2=1)[:,:,np.newaxis] #(ngllx,1,ngllz) * (ngllx,ngllx,1) = (ngllx,ngllx,ngllz)
        dux_dxi=sum(temp_dux_dxi) #du_x/dxi, (ngllx,ngllz)

        temp_dux_dxi=displ[1,ibool[:,:,i]][:,np.newaxis]*np.swapaxes(a=hprime_xx,axis1=0,axis2=1)[:,:,np.newaxis]
        duz_dxi=sum(temp_dux_dxi) #du_z/dxi

        temp_dux_dgamma=displ[0,ibool[:,:,i]][:,:,np.newaxis]*np.swapaxes(a=hprime_zz,axis1=0,axis2=1) #(ngllx,ngllz,1) * (ngllz,ngllz) = (ngllx,ngllz,ngllz)
        dux_dgamma=np.sum(a=temp_dux_dgamma,axis=1) #du_x/dgamma, (ngllx,ngllz)

        temp_dux_dgamma=displ[1,ibool[:,:,i]][:,:,np.newaxis]*np.swapaxes(a=hprime_zz,axis1=0,axis2=1)
        duz_dgamma=np.sum(a=temp_dux_dgamma,axis=1) #du_z/dgamma

        dux_dxl=dux_dxi*xix[:,:,i]+dux_dgamma*gammax[:,:,i] #du_x/dx, local
        dux_dzl=dux_dxi*xiz[:,:,i]+dux_dgamma*gammaz[:,:,i] #du_x/dz, local

        duz_dxl=duz_dxi*xix[:,:,i]+duz_dgamma*gammax[:,:,i] #du_z/dx, local
        duz_dzl=duz_dxi*xiz[:,:,i]+duz_dgamma*gammaz[:,:,i] #du_z/dz, local

        sigma_xx=lambdaplus2mu*dux_dxl+lambda_*duz_dzl
        sigma_xz=mu*(duz_dxl+dux_dzl)
        sigma_zz=lambdaplus2mu*duz_dzl+lambda_*dux_dxl
        sigma_zx=sigma_xz

        tempx1=wzgll*jacobian[:,:,i]*(sigma_xx*xix[:,:,i]+sigma_zx*xiz[:,:,i]) #(ngllz) * (ngllx,ngllz) = (ngllx,ngllz)
        tempz1=wzgll*jacobian[:,:,i]*(sigma_xz*xix[:,:,i]+sigma_zz*xiz[:,:,i])

        tempx2=wxgll[:,np.newaxis]*jacobian[:,:,i]*(sigma_xx*gammax[:,:,i]+sigma_zx*gammaz[:,:,i]) #(ngllx,1) * (ngllx,ngllz) = (ngllx,ngllz)
        tempz2=wxgll[:,np.newaxis]*jacobian[:,:,i]*(sigma_xz*gammax[:,:,i]+sigma_zz*gammaz[:,:,i])
        
        tempx3=tempx1[:,np.newaxis]*hprimewgll_xx[:,:,np.newaxis] #(ngllx,1,ngllz) * (ngllx,ngllx,1) = (ngllx,ngllx,ngllz)
        tempx4=np.swapaxes(a=tempx2,axis1=0,axis2=1)[:,:,np.newaxis]*hprimewgll_zz[:,np.newaxis] #(ngllz,ngllx,1) * (ngllz,1,ngllz) = (ngllz,ngllx,ngllz)
        accel[0,ibool[:,:,i]]-=sum(tempx3)+sum(tempx4) #(ngllx,ngllz) + (ngllx,ngllz) = (ngllx,ngllz)

        tempz3=tempz1[:,np.newaxis]*hprimewgll_xx[:,:,np.newaxis]
        tempz4=np.swapaxes(a=tempz2,axis1=0,axis2=1)[:,:,np.newaxis]*hprimewgll_zz[:,np.newaxis]
        accel[1,ibool[:,:,i]]-=sum(tempz3)+sum(tempz4)

    accel[1,ibool[igll_source,jgll_source,nspec_source]]-=factor_amplitude*(1-2*a*(it*deltat-t0)**2)*np.exp(-a*(it*deltat-t0)**2)
    accel*=rmass_inverse #(ndim,nglob) * (nglob) = (ndim,nglob)
    veloc+=newmark_gamma*deltat*accel

    seismogram.append(displ[1,ibool[igll_receiver,jgll_receiver,nspec_receiver]])

with open(file="seismogram.txt",mode="w") as file:
    for it in range(nstep):
        file.write(f"{it*deltat} {seismogram[it]}\n")
