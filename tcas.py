from obspy import read

import numpy as np

nsi=10 #no. of stacking iterations
pjgl=1.2 #stack index (the Lp norm used)
erl=1.0095 #error factor
emin=0.01 #minimum error
emax=0.15 #maximum error
stkwb=34 #start of stack window
stkwl=13 #length of stack window
dtcmin=-0.5 #min bound on differential time search
dtcmax=0.5 #max bound on differential time search

#user0 header=moveout correction by ak135
#npts don't have to be same, but begin time must all be same
st=read("*.sac.ft.cut")
delta=st[0].stats.delta

dtmo=np.array([-st[i].stats.sac.user0 for i in range(len(st))]) #moveout correction
imo=np.rint(dtmo/delta).astype(int)
dtcs=np.zeros(len(st)) #local time shift from stacking iteration
nstkwb=int(stkwb/delta) #no. of samples to start of stack window
nstkwl=int(stkwl/delta) #no. of samples in stack window
jim1=int(np.rint(dtcmin/delta))
jim2=int(np.rint(dtcmax/delta))
errs=np.empty(len(st)) #pick error estimated from trace power

#normalization
for i in range(len(st)):
    st[i].data/=max(abs(st[i].data))

#initial stacking
zssl=np.zeros(nstkwl+1) #linear trace stack
zscp=np.zeros(nstkwl+1) #quadratic trace stack
lmn=np.rint((-dtmo+dtcs)/delta).astype(int)

for i in range(len(st)):
    zssl+=st[i].data[nstkwb-lmn[i]-1:nstkwb+nstkwl-lmn[i]]
    zscp+=st[i].data[nstkwb-lmn[i]-1:nstkwb+nstkwl-lmn[i]]**2

pstakn=sum(zscp)/(len(st)*stkwl) #L2 measure of trace misfit
print(f"initial pstakn={pstakn}")
zssl/=max(abs(zssl))

#adaptive stacking
for i in range(nsi):
    for j in range(len(st)):
        wsp=np.empty(jim2-jim1+1) #power of weighted stack

        for k in range(jim2-jim1+1):
            wsp[k]=sum(abs(zssl-st[j].data[nstkwb+imo[j]+(k+jim1)-1:nstkwb+nstkwl+imo[j]+(k+jim1)])**pjgl)/stkwl #Lp norm misfit

        jm=np.argmin(wsp) #minimun value's index

        #error estimation
        if i==nsi-1: #only last iteration
            wm=np.min(wsp) #minimum value
            swl,swr=0,0 #switches for half width error determination

            if 0<jm: #minimum index is not left boundry
                for l in range(jm-1,-1,-1):
                    if wm*erl<wsp[l]: #found intersection in left range
                        swl=1
                        errl=((jm-l)-(wsp[l]-wm*erl)/(wsp[l]-wsp[l+1]))*delta #linear interpolation
                        break

            if jm<jim2-jim1: #minimum index is not right boundary
                for l in range(jm+1,jim2-jim1):
                    if wm*erl<wsp[l]: #found intersection in right range
                        swr=1
                        errr=((l-jm)-(wsp[l]-wm*erl)/(wsp[l]-wsp[l-1]))*delta #linear interpolation
                        break

            if swl==1 and swr==1: #left, right error exist
                err=0.5*(errr+errl) #average
            elif swl==1: #only left error exists
                err=errl
            elif swr==1: #only right error exists
                err=errr
            else: #neither
                err=emax

            if err<emin: #bound to minimum error
                err=emin
            elif emax<err: #bound to maximum error
                err=emax

            errs[j]=err

        dtcs[j]=-(jm+jim1)*delta #new time shift

    #stacking
    zssl=np.zeros(nstkwl+1)
    zscp=np.zeros(nstkwl+1)
    lmn=np.rint((-dtmo+dtcs)/delta).astype(int)

    for j in range(len(st)):
        zssl+=st[j].data[nstkwb-lmn[j]-1:nstkwb+nstkwl-lmn[j]]
        zscp+=st[j].data[nstkwb-lmn[j]-1:nstkwb+nstkwl-lmn[j]]**2

    pstakn=sum(zscp)/(len(st)*stkwl)
    print(f"{i+1} iteration pstakn={pstakn}")
    zssl/=max(abs(zssl))
