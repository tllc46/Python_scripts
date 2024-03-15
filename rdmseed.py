from datetime import datetime,timedelta
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt

def nibble_split(binary):
    nibbles=[]
    nibbles.append(binary&3)
    for i in range(15):
        binary>>=2
        nibbles.append(binary&3)
    nibbles.reverse()
    return nibbles

def data_split(binary,size):
    num=30//size
    differences=[]
    difference=binary&((1<<size)-1)
    if difference>>(size-1):
        difference-=1<<size
    differences.append(difference)
    for i in range(num-1):
        binary>>=size
        difference=binary&((1<<size)-1)
        if difference>>(size-1):
            difference-=1<<size
        differences.append(difference)
    differences.reverse()
    return differences

def decode_data(file):
    difference=[]
    global data
    break_i=False
    for i in range(frame_count):
        nibble=nibble_split(int.from_bytes(file.read(4)))
        for j in range(1,16): #j=0 is for nibbles
            if nibble[j]==0:
                if i==0 and j==1: #first data
                    x0=int.from_bytes(file.read(4),signed=True)
                elif i==0 and j==2: #last data
                    xn=int.from_bytes(file.read(4),signed=True)
                else: #end of data
                    file.seek((16-j)*4+(frame_count-1-i)*64,1) #move to end of record
                    break_i=True
                    break
            elif nibble[j]==1:
                difference+=list(unpack(">4b",file.read(4)))
            elif nibble[j]==2:
                binary=int.from_bytes(file.read(4))
                decode_nibble=binary>>30
                if decode_nibble==1:
                    difference+=data_split(binary,30)
                elif decode_nibble==2:
                    difference+=data_split(binary,15)
                else: #decode_nibble==3:
                    difference+=data_split(binary,10)
            else: #nibble[j]==3:
                binary=int.from_bytes(file.read(4))
                decode_nibble=binary>>30
                if decode_nibble==0:
                    difference+=data_split(binary,6)
                elif decode_nibble==1:
                    difference+=data_split(binary,5)
                else: #decode_nibble==2:
                    difference+=data_split(binary,4)
        if break_i:
            break
    difference=np.array(difference[1:])
    data=np.concatenate(([x0],x0+np.cumsum(difference)))
    return None

activity_flags=["calibration signals present",
        "time correction applied",
        "beginning of an event, station triggers",
        "end of the event, station detriggers",
        "+ leap second included",
        "- leap second included",
        "event in progress"]
IO_and_clock_flags=["station volume parity error possibly present",
        "long record read (possibly no problem)",
        "short record read (record padded)",
        "start of time series",
        "end of time series",
        "clock locked"]
data_quality_flags=["amplifier saturation detected",
        "digitizer clipping detected",
        "spikes detected",
        "glitches detected",
        "missing/padded data present",
        "telemetry synchronization error",
        "a digital filter may be charging",
        "time tag is questionable"]
encoding_format={0:"text, UTF-8 allowed, use ASCII for maximum portability, no structure defined",
        1:"16 bit integer (2’s complement)",
        3:"32 bit integer (2’s complement)",
        4:"32 bit floating point (IEEE float)",
        5:"64 bit double precesion floating point (IEEE double)",
        10:"Steim-1 compression",
        11:"Steim-2 compression",
        19:"Steim-3 compression",
        100:"opaque data, only for use in special scenarios, not intended for archiving"}
blockette_type={0:"none",
        1000:"data only SEED (1000)",
        1001:"data extension (1001)"}
word_order=["little endian",
        "big endian"]

with open("AU.WB9..BHZ",mode="rb") as file:
    traces=[]
    times=[]
    #first record
    #read and parse full header
    header_word=file.read(20).decode("ascii")
    header_num=unpack(">HHBBBBH HHHBBBBlHH",file.read(28))
    zero_date=datetime.strptime(f"{header_num[0]} {header_num[1]} {header_num[2]} {header_num[3]} {header_num[4]} {header_num[6]*100:06}","%Y %j %H %M %S %f") #starting date
    npts=header_num[7] #number of samples
    delta=1/header_num[8] #delta

    b1000=unpack(">HHBBBB",file.read(8)) #blockette 1000
    frame_count=(1<<(b1000[4]-6))-1 #frame count

    file.seek(8,1) #dummy bytes

    #decode data
    decode_data(file)
    time=np.arange(npts)*delta
    traces.append(data)
    times.append(time)
    end_date=zero_date+timedelta(seconds=(npts-1)*delta)

    while file.read(1):
        #read header
        file.seek(19,1)
        sub_header=unpack(">HHBBBBH H",file.read(12)) #parse only starting date and number of samples
        start_date=datetime.strptime(f"{sub_header[0]} {sub_header[1]} {sub_header[2]} {sub_header[3]} {sub_header[4]} {sub_header[6]*100:06}","%Y %j %H %M %S %f") #starting date
        npts=sub_header[7] #number of samples
        file.seek(32,1)

        #decode data
        decode_data(file)
        offset=(start_date-zero_date).total_seconds()
        time=offset+np.arange(npts)*delta
        print(len(time),npts,delta)
        if delta<(start_date-end_date).total_seconds(): #gap
            traces.append(data)
            times.append(time)
        else: #no gap
            traces[-1]=np.concatenate((traces[-1],data))
            times[-1]=np.concatenate((times[-1],time))
        end_date=start_date+timedelta(seconds=(npts-1)*delta)

#print header information
print(f"sequence number={header_word[:6]}")
print(f"data header/quality indicator={header_word[6]}")
print(f"reserved byte={header_word[7]}")
print(f"station identifier code={header_word[8:13]}")
print(f"location identifier={header_word[13:15]}")
print(f"channel identifier={header_word[15:18]}")
print(f"network code={header_word[18:]}")
print(f"year={header_num[0]}")
print(f"day of year={header_num[1]}")
print(f"hour={header_num[2]}")
print(f"minute={header_num[3]}")
print(f"second={header_num[4]}")
print(f"unused={header_num[5]}")
print(f"second(10000)={header_num[6]}")
print(f"number of samples={header_num[7]}")
print(f"sample rate factor={header_num[8]}")
print(f"sample rate multiplier={header_num[9]}")
print("activity flags=",activity_flags[header_num[10]])
print("I/O and clock flags=",IO_and_clock_flags[header_num[11]])
print("data quality flags=",data_quality_flags[header_num[12]])
print(f"number of blockettes that follow={header_num[13]}")
print(f"time correction={header_num[14]}")
print(f"beginning of data={header_num[15]}")
print(f"fisrt blockette={header_num[16]}")
print() #blockette 1000
print("blockette type=",blockette_type[b1000[0]])
print(f"next blockette's byte number={b1000[1]}")
print("encoding format=",encoding_format[b1000[2]])
print("word order=",word_order[b1000[3]])
print(f"data record length={b1000[4]}")
print(f"reserved byte={b1000[5]}")

#plot waveform
fig=plt.figure()
ax=fig.subplots()
for i in range(len(traces)):
    ax.plot(times[i],traces[i],color="black")
fig.savefig("trace.png")
