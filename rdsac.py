from struct import unpack

import numpy as np
import matplotlib.pyplot as plt

enum=[0,
        "time series", #file type
        "real, imaginary spectral", #file type
        "amplitude, phase spectral", #file type
        "general x vs. y", #file type
        "unknown 1", #file type
        "displacement(nm)", #dependent variable type
        "velocity(nm/sec)", #dependent variable type
        "acceleration(nm/sec^2)", #dependent variable type
        "begin time", #zero time type
        "midnight of zero time day", #zero time type
        "origin time", #zero time type
        "arrival time", #zero time type
        "time pick 1", #zero time type
        "time pick 2", #zero time type
        "time pick 3", #zero time type
        "time pick 4", #zero time type
        "time pick 5", #zero time type
        "time pick 6", #zero time type
        "time pick 7", #zero time type
        "time pick 8", #zero time type
        "time pick 9", #zero time type
        "time pick 10", #zero time type
        "radial(Nevada test site)", #component type
        "tangential(Nevada test site)", #component type
        "radial(event)", #component type
        "tangential(event)", #component type
        "north positive", #component type
        "east positive", #component type
        "horizontal(arbitrary)", #component type
        "down positive", #component type
        "up positive", #component type
        "Lawrence Livermore laboratory broadband", #instrument type
        "world wide standardized seismograph network 15-100", #instrument type
        "world wide standardized seismograph network 30-100", #instrument type
        "high gain, long period", #instrument type
        "seismic research observatory", #instrument type
        "nuclear event 1", #event type
        "nuclear pre-shot event", #event type
        "nuclear post-shot event", #event type
        "earthquake 1", #event type
        "foreshock", #event type
        "aftershock", #event type
        "chemical explosion", #event type
        "other 1", #event type
        "good", #data quality
        "glitches", #data quality
        "dropouts", #data quality
        "low signal to noise ratio", #data quality
        "real data", #synthetic data flag
        "velocity(V)", #dependent variable type
        "general xy vs. z", #file type
        "body wave magnitude", #magnitude type
        "surface wave magnitude", #magnitude type
        "local magnitude", #magnitude type
        "moment magnitude", #magnitude type
        "duration magnitude", #magnitude type
        "user defined magnitude", #magnitude type
        "national earthquake information center", #magnitude information source
        "quick preliminary determination of epicenter", #magnitude information source
        "weekly preliminary determination of epicenter", #magnitude information source
        "preliminary determination of epicenter", #magnitude information source
        "internation seismological center", #magnitude information source
        "reviewed event bulletin", #magnitude information source
        "US geological survey", #magnitude information source
        "UC Berkeley", #magnitude information source
        "California institute of technology", #magnitude information source
        "Lawrence Livermore national laboratory", #magnitude information source
        "event location(computer program)", #magnitude information source
        "joint seismic observation program", #magnitude information source
        "user", #magnitude information source
        "unknown 2", #magnitude information source
        "quarry/mine blast confirmed by quarry", #event type
        "quarry/mine blast with designed shot info-ripple fired", #event type
        "quarry/mine blast with observed shot info-ripple fired", #event type
        "quarry/mine blast: single shot", #event type
        "quarry/mining-induced events: tremors, rockbursts", #event type
        "earthquake 2", #event type
        "earthquakes in a swarm/aftershock sequence", #event type
        "felt earthquake", #event type
        "marine explosion", #event type
        "other explosion", #event type
        "nuclear event 2", #event type
        "nuclear cavity collapse", #event type
        "other 2", #event type
        "unknown(local event)", #event type
        "unknown(regional event)", #event type
        "unknown(teleseismic event)", #event type
        "undetermined/conflicting information", #event type
        "damaging earthquake", #event type
        "probable earthquake", #event type
        "probable explosion", #event type
        "mine collapse", #event type
        "probable mine blast", #event type
        "geyser", #event type
        "lightning", #event type
        "meteroic event", #event type
        "odors", #event type
        "Sun", #body type
        "Mercury", #body type
        "Venus", #body type
        "Earth", #body type
        "Moon", #body type
        "Mars"] #body type

logic=[False,True]

file=open(file="NS.N068..HHE.D.2021.280.000000.sac.swap",mode="rb")
fhdr=unpack(">70f",file.read(70*4))
delta=fhdr[0]
nhdr=unpack(">40i",file.read(40*4))
npts=nhdr[9]
khdr=file.read(192).decode(encoding="ascii")
data=unpack(f">{npts}f",file.read(npts*4))
file.close()

float_header=list(fhdr)
int_header=list(nhdr)
data=list(data)
time=np.arange(stop=npts)*delta

#실수 header
for i in range(70):
    if float_header[i]==-12345:
        float_header[i]="undefined"

#정수 header
for i in range(40):
    if int_header[i]==-12345:
        int_header[i]="undefined"
    else:
        if 15<=i<=34: #열거 header
            int_header[i]=enum[int_header[i]]
        elif 35<=i: #논리 header
            int_header[i]=logic[int_header[i]]

#문자 header
char_header=[]
for i in range(23):
    if i==0:
        if "\0" in khdr[:8]:
            char_header.append(khdr[:khdr[:8].index("\0")])
        else:
            char_header.append(khdr[:8])
    elif i==1:
        if "\0" in khdr[8:24]:
            char_header.append(khdr[8:8+khdr[8:24].index("\0")])
        else:
            char_header.append(khdr[8:24])
    else:
        if "\0" in khdr[8+8*i:16+8*i]:
            char_header.append(khdr[8+8*i:8+8*i+khdr[8+8*i:16+8*i].index("\0")])
        else:
            char_header.append(khdr[8+8*i:16+8*i])
    if char_header[i][:6]=="-12345":
        char_header[i]="undefined"

#header 정보 출력
print("<required fields>")
print(f"number of points={int_header[9]}")
print(f"header version={int_header[6]}")
print(f"beginning value={float_header[5]}")
print(f"ending value={float_header[6]}")
print(f"type of file={int_header[15]}")
print(f"evenly spaced={int_header[35]}")
print(f"delta={float_header[0]}")
print()
print("<time fields>")
print(f"observed delta={float_header[4]}")
print(f"type of dependent variable={int_header[16]}")
print(f"minimum of dependent variable={float_header[1]}")
print(f"maximum of dependent variable={float_header[2]}")
print(f"mean of dependent variable={float_header[56]}")
print(f"zero time year={int_header[0]}")
print(f"zero time juldian day={int_header[1]}")
print(f"zero time hour={int_header[2]}")
print(f"zero time minute={int_header[3]}")
print(f"zero time second={int_header[4]}")
print(f"zero time millisecond={int_header[5]}")
print(f"type of zero time={int_header[17]}")
print(f"origin time={float_header[7]}")
print("origin time identification=",char_header[3])
print(f"stored number of points={int_header[10]}")
print(f"stored begin={float_header[54]}")
print(f"stored delta={float_header[55]}")
print()
print("<phase picks>")
print(f"arrival time={float_header[8]}")
print("arrival time identification=",char_header[4])
print(f"finish time={float_header[20]}")
print("finish time identification=",char_header[15])
print(f"time pick 1={float_header[10]}")
print(f"time pick 2={float_header[11]}")
print(f"time pick 3={float_header[12]}")
print(f"time pick 4={float_header[13]}")
print(f"time pick 5={float_header[14]}")
print(f"time pick 6={float_header[15]}")
print(f"time pick 7={float_header[16]}")
print(f"time pick 8={float_header[17]}")
print(f"time pick 9={float_header[18]}")
print(f"time pick 10={float_header[19]}")
print("time pick 1 identification=",char_header[5])
print("time pick 2 identification=",char_header[6])
print("time pick 3 identification=",char_header[7])
print("time pick 4 identification=",char_header[8])
print("time pick 5 identification=",char_header[9])
print("time pick 6 identification=",char_header[10])
print("time pick 7 identification=",char_header[11])
print("time pick 8 identification=",char_header[12])
print("time pick 9 identification=",char_header[13])
print("time pick 10 identification=",char_header[14])
print()
print("<instrument fields>")
print("name of instrument=",char_header[22])
print(f"type of instrument={int_header[19]}")
print(f"instrument response 1={float_header[21]}")
print(f"instrument response 2={float_header[22]}")
print(f"instrument response 3={float_header[23]}")
print(f"instrument response 4={float_header[24]}")
print(f"instrument response 5={float_header[25]}")
print(f"instrument response 6={float_header[26]}")
print(f"instrument response 7={float_header[27]}")
print(f"instrument response 8={float_header[28]}")
print(f"instrument response 9={float_header[29]}")
print(f"instrument response 10={float_header[30]}")
print()
print("<station fields>")
print("network=",char_header[20])
print("station name=",char_header[0])
print(f"station region={int_header[20]}")
print(f"station latitude={float_header[31]}")
print(f"station longitude={float_header[32]}")
print(f"station elevation={float_header[33]}")
print(f"station depth={float_header[34]}")
print(f"component azimuth={float_header[57]}")
print(f"component incident angle={float_header[58]}")
print("component name=",char_header[19])
print(f"positive polarity={int_header[36]}")
print()
print("<event fields>")
print("event name=",char_header[1])
print(f"event region={int_header[21]}")
print(f"event latitude={float_header[35]}")
print(f"event longitude={float_header[36]}")
print(f"event elevation={float_header[37]}")
print(f"event depth={float_header[38]}")
print(f"magnitude={float_header[39]}")
print(f"type of magnitude={int_header[25]}")
print(f"source of magnitude information={int_header[26]}")
print(f"type of event={int_header[22]}")
print(f"event ID={int_header[8]}")
print(f"origin ID={int_header[7]}")
print(f"waveform ID={int_header[11]}")
print("location identifier=",char_header[2])
print(f"distance={float_header[50]}")
print(f"azimuth={float_header[51]}")
print(f"back azimuth={float_header[52]}")
print(f"great circle arc length={float_header[53]}")
print(f"body definition={int_header[27]}")
print()
print("<miscellaneous fields>")
print(f"calculate distance, azimuth={int_header[38]}")
print(f"quality of data={int_header[23]}")
print(f"synthetic data flag={int_header[24]}")
print("date data read=",char_header[21])
print(f"user defined variable 1={float_header[40]}")
print(f"user defined variable 2={float_header[41]}")
print(f"user defined variable 3={float_header[42]}")
print(f"user defined variable 4={float_header[43]}")
print(f"user defined variable 5={float_header[44]}")
print(f"user defined variable 6={float_header[45]}")
print(f"user defined variable 7={float_header[46]}")
print(f"user defined variable 8={float_header[47]}")
print(f"user defined variable 9={float_header[48]}")
print(f"user defined variable 10={float_header[49]}")
print("user defined identification 1=",char_header[16])
print("user defined identification 2=",char_header[17])
print("user defined identification 3=",char_header[18])
print(f"okay to overwrite={int_header[37]}")
print(f"size of x={int_header[12]}")
print(f"size of y={int_header[13]}")
print(f"minimum of x={float_header[59]}")
print(f"maximum of x={float_header[60]}")
print(f"minimum of y={float_header[61]}")
print(f"maximum of y={float_header[62]}")
print()
print("<unused fields>")
print(f"multiplying scale factor={float_header[3]}")
print(f"format={float_header[9]}")
print(f"time adjustment={float_header[63]}")
print(f"floating header 65={float_header[64]}")
print(f"floating header 66={float_header[65]}")
print(f"floating header 67={float_header[66]}")
print(f"floating header 68={float_header[67]}")
print(f"floating header 69={float_header[68]}")
print(f"floating header 70={float_header[69]}")
print(f"integer header 15={int_header[14]}")
print(f"component={int_header[18]}")
print(f"enumerated header 14={int_header[28]}")
print(f"enumerated header 15={int_header[29]}")
print(f"enumerated header 16={int_header[30]}")
print(f"enumerated header 17={int_header[31]}")
print(f"enumerated header 18={int_header[32]}")
print(f"enumerated header 19={int_header[33]}")
print(f"enumerated header 20={int_header[34]}")
print(f"logical header 5={int_header[39]}")

#파형 그림
fig=plt.figure()
ax=fig.subplots()
ax.plot(time,data,color="black")
fig.savefig(fname="trace.png")
