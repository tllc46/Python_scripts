#astack/tcas.f 참고

from obspy import read

import numpy as np

nsi=10 #stack 반복 횟수
pjgl=1.2 #Lp norm 지수
erl=1.0095 #오차 factor
emin=0.01 #최소 오차
emax=0.15 #최대 오차
stkwb=34 #stack 구간 시작 시각
stkwl=13 #stack 구간 길이
dtcmin=-0.5 #최소 time shift
dtcmax=0.5 #최대 time shift

#user0 header는 ak135 모형에 근거한 초기 moveout 저장
#상대적으로 나중에 도착하면 양수, 먼저 도착하면 음수
#npts는 같을 필요는 없지만, 시작 시각은 모두 같아야 한다
st=read(pathname_or_url="*.sac.ft.cut",format="SAC",byteorder="little")
delta=st[0].stats.delta

dtmo=np.array(object=[-st[i].stats.sac.user0 for i in range(len(st))]) #ak135 moveout
imo=int(np.rint(dtmo/delta))
dtcs=np.zeros(shape=len(st)) #1회 stack으로 결정된 time shift
#ak135에 대해 나중에 도착하면 음수, 먼저 도착하면 양수
nstkwb=int(stkwb/delta) #data 시작부터 stack 구간 시작까지 npts
nstkwl=int(stkwl/delta) #stack 구간 npts
jim1=int(np.rint(dtcmin/delta))
jim2=int(np.rint(dtcmax/delta))
errs=np.empty(shape=len(st)) #오차

#정규화
for i in range(len(st)):
    st[i].data/=max(abs(st[i].data))

#초기 stack
zssl=np.zeros(shape=nstkwl+1) #선형 stack
zscp=np.zeros(shape=nstkwl+1) #이차 stack
lmn=np.rint((-dtmo+dtcs)/delta).astype(dtype=int)

for i in range(len(st)):
    zssl+=st[i].data[nstkwb-lmn[i]-1:nstkwb+nstkwl-lmn[i]]
    zscp+=st[i].data[nstkwb-lmn[i]-1:nstkwb+nstkwl-lmn[i]]**2

pstakn=sum(zscp)/(len(st)*stkwl)
print(f"initial pstakn={pstakn}")
zssl/=max(abs(zssl))

#adaptive stacking
for i in range(nsi):
    for j in range(len(st)):
        wsp=np.empty(shape=jim2-jim1+1) #time shift 탐색 구간에서 Lp norm

        for k in range(jim2-jim1+1):
            wsp[k]=sum(abs(zssl-st[j].data[nstkwb+imo[j]+(k+jim1)-1:nstkwb+nstkwl+imo[j]+(k+jim1)])**pjgl)/stkwl #Lp norm

        jm=np.argmin(a=wsp) #Lp norm이 최소가 되는 색인

        #오차
        if i==nsi-1: #마지막 반복 단계에서만 계산
            wm=wsp[jm] #최소 Lp norm
            swl,swr=0,0 #음/양의 오차 존재 여부

            if 0<jm: #최소 Lp norm이 탐색 구간 시작 이후 존재
                for l in range(jm-1,-1,-1): #음의 오차 탐색
                    if wm*erl<wsp[l]: #오차 교차점 발견
                        swl=1 #음의 오차 존재
                        errl=((jm-l)-(wsp[l]-wm*erl)/(wsp[l]-wsp[l+1]))*delta #선형 보간
                        break

            if jm<jim2-jim1: #최소 Lp norm이 탐색 구간 끝 이전 존재
                for l in range(jm+1,jim2-jim1): #양의 오차 탐색
                    if wm*erl<wsp[l]: #오차 교차점 발견
                        swr=1 #양의 오차 존재
                        errr=((l-jm)-(wsp[l]-wm*erl)/(wsp[l]-wsp[l-1]))*delta #선형 보간
                        break

            if swl==1 and swr==1: #음/양의 오차 모두 존재
                err=0.5*(errr+errl) #평균
            elif swl==1: #음의 오차만 존재
                err=errl
            elif swr==1: #양의 오차만 존재
                err=errr
            else: #아무 오차도 없다
                err=emax

            if err<emin: #최소 오차보다 작은 경우
                err=emin
            elif emax<err: #최대 오차보다 큰 경우
                err=emax

            errs[j]=err

        dtcs[j]=-(jm+jim1)*delta #time shift 결정

    #stacking
    zssl=np.zeros(shape=nstkwl+1)
    zscp=np.zeros(shape=nstkwl+1)
    lmn=np.rint((-dtmo+dtcs)/delta).astype(dtype=int)

    for j in range(len(st)):
        zssl+=st[j].data[nstkwb-lmn[j]-1:nstkwb+nstkwl-lmn[j]]
        zscp+=st[j].data[nstkwb-lmn[j]-1:nstkwb+nstkwl-lmn[j]]**2

    pstakn=sum(zscp)/(len(st)*stkwl)
    print(f"{i+1} iteration pstakn={pstakn}")
    zssl/=max(abs(zssl))
