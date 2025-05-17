#sac-102.0/utils/sgftops.c 참고

from os import getcwd,geteuid
from os.path import getmtime
from pwd import getpwuid
from time import ctime
from struct import unpack

import numpy as np
import pandas as pd

scale_option=False
id_option=False
width_in=1
colortable=False

def EndPath():
    global fillFlag

    if not fillFlag:
        ofp.write("s\n")
    else:
        ofp.write("fill\n")
        ofp.write("0 g\n")
        fillFlag=True

def execute_buffer():
    global textangle,imagecount,icount,doing_image
    global needmove,fillFlag
    global done

    #are we in the middle of writing out an image?
    if doing_image:
        #print out rgb values 24 to a line
        bufcount=0
        hexnums=cbuf2.hex()

        while icount+24<=imagecount and bufcount+12<=buflen:
            ofp.write(hexnums[4*bufcount:4*(bufcount+12)]+"\n")
            icount+=24
            bufcount+=12

        if bufcount<buflen:
            nwrite=min(2*(buflen-bufcount),imagecount-icount)
            ofp.write(hexnums[4*bufcount:4*bufcount+2*nwrite]+"\n")
            icount+=nwrite

        if icount==imagecount:
            doing_image=False
            ofp.write(" grestore\n")

        return

    i=0
    while i<buflen and not done:
        if needmove and (0<=buffer[i] or buffer[i]==-5 or buffer[i]==-10):
            EndPath()
            ofp.write(f"{x} {y} m\n")
            needmove=False

        if 0<=buffer[i]:
            x=buffer[i]
            y=buffer[i+1]
            ofp.write(f"{x} {y} L\n")
            i+=2

        elif buffer[i]==-1: #No-op command
            i+=1

        elif buffer[i]==-2: #End of plot command
            done=True
            i+=2

        elif buffer[i]==-3: #Move command
            x=buffer[i+2]
            y=buffer[i+3]
            needmove=True
            i+=4

        elif buffer[i]==-4: #Color command - default SAC color map
            indx=buffer[i+1]
            if indx==1:
                newindx=buffer[i+2]
                if 0<=newindx<ict and newindx!=indx:
                    indx=newindx
                    EndPath()
                    ofp.write(f"{ctab[indx,0]} {ctab[indx,1]} {ctab[indx,2]} setrgbcolor\n")
                else:
                    EndPath()
                    ofp.write("0 0 0 setrgbcolor\n")
                    print("Illegal color value encountered - set to default color black")
                i+=3
            else:
                EndPath()
                ofp.write(f"{buffer[i+2]/255} ")
                ofp.write(f"{buffer[i+3]/255} ")
                ofp.write(f"{buffer[i+4]/255} setrgbcolor\n")
                i+=5

        elif buffer[i]==-5: #Hardware text command
            y=buffer[i+7]
            string=cbuf2[2*(i+8):2*(i+8)+y].decode()
            ofp.write("gsave\n")
            if textangle:
                ofp.write(f"{textangle} rotate\n")
                ofp.write("32000 10 72 mul div 24000 7.5 72 mul div scale\n")
                ofp.write("({string}) show\n")
                ofp.write("grestore\n")
            i+=8+(y+1)//2

        elif buffer[i]==-6: #Hardware text size
            i+=4 #(ignored)

        elif buffer[i]==-7: #Linestyle command
            EndPath()
            if buffer[i+2]<=0:
                buffer[i+2]=1
            ofp.write(f"{linemodes[(buffer[i+2]-1)%numlstyls]}\n")
            i+=3

        elif buffer[i]==-9: #Linewidth command
            EndPath()
            if buffer[i+2]<=0:
                buffer[i+2]=1
            new_width=buffer[i+2]
            if new_width<1 or 10<new_width:
                new_width=1
            new_width*=base_width
            ofp.write(f"{new_width} setlinewidth\n")
            i+=3

        elif buffer[i]==-10: #polyfill command
            igray=buffer[i+2]
            if not igray:
                gray_value=1
            elif igray==1:
                gray_value=0.8
            else:
                gray_value=0
            #move command
            x=buffer[i+5]
            y=buffer[i+6]
            EndPath()
            ofp.write(f"{gray_value} g\n")
            ofp.write(f"{x} {y} m\n")
            fillFlag=True
            i+=7

        elif buffer[i]==-11: #plot filled rectangle
            x=buffer[i+2]
            y=buffer[i+3]
            dx=buffer[i+4]
            dy=buffer[i+5]
            indx=buffer[i+6]
            if indx<0 or ict<=indx:
                print("Illegal color value: set to last index")
                indx=ict-1
            ofp.write(f"{x} {y} {dx} {dy} {ctab[indx,0]} {ctab[indx,1]} {ctab[indx,2]} F\n")
            i+=7

        elif buffer[i]==-12: #set text angle
            textangle=buffer[i+2]
            i+=3

        elif buffer[i]==-13: #color image
            iwidth=buffer[i+1]
            iheight=buffer[i+2]
            xloc=buffer[i+3]
            yloc=buffer[i+4]
            ofp.write(" gsave\n")
            ofp.write(f" /picstr {3*iwidth} string def\n")
            ofp.write(f" {xloc} {yloc} translate\n")
            ofp.write(f" {32*iwidth} {32*iheight} scale\n")
            ofp.write(f" {iwidth} {iheight} 8\n")
            ofp.write(f" [ {iwidth} 0 0 {-iheight} 0 {iheight} ]\n")
            ofp.write(" {currentfile picstr readhexstring pop}\n")
            ofp.write(" false 3\n")
            ofp.write(" colorimage\n\n")
            imagecount=iwidth*iheight*3
            doing_image=True
            icount=0
            i+=5

        elif buffer[i]==-14: #RGB fill
            ofp.write("fill\n")
            i+=1

        else:
            count=buffer[i+1]
            print(f"unknown op code {buffer[i]}, count = {count}")
            i+=2+count

numlstyls=10
linemodes=["ls1","ls2","ls3","ls4","ls5","ls6","ls7","ls8","ls9","ls10"]
textangle=0
doing_image=False

needmove=False
fillFlag=False

#main
done=False

#open the sgf file
file_in=open(file="f001.sgf",mode="rb")
#open the postscriptfile to write
ofp=open(file="f003.ps",mode="w")

#get the info for the plot id label
usnam=getpwuid(geteuid()).pw_name
path=getcwd()
dattime=ctime(getmtime(filename="f001.sgf"))

#use optional line width if argument present
base_width=int(40*width_in)

if colortable:
    #user defined colortable
    df=pd.read_csv(filepath_or_buffer=colorfile,sep=" ",names=["red","green","blue"])
    ctab=df.to_numpy()
else:
    #default colortable
    ctab=np.array(object=[[1,1,1],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,0,1],[0,0,0]])
ict=len(ctab)

#define some abreviations for some postscript commands
ofp.write("%!PS-Adobe-3.0 EPSF-3.0\n")
ofp.write("%%BoundingBox: 0 0 612 792\n")
ofp.write("%%Creator: sgftops (20090313)\n")
ofp.write("%%Pages: 1\n")
ofp.write("%%EndComments\n")
ofp.write("%%BeginProlog\n")
ofp.write("/Helvetica findfont 11 scalefont setfont\n")
ofp.write("/ls1 {[] 0 setdash} def\n")
ofp.write("/ls2 {[200 150] 0 setdash} def\n")
ofp.write("/ls3 {[600 150] 0 setdash} def\n")
ofp.write("/ls4 {[600 150 75 150] 0 setdash} def\n")
ofp.write("/ls5 {[600 150 200 150] 0 setdash} def\n")
ofp.write("/ls6 {[600 150 75 150 75 150] 0 setdash} def\n")
ofp.write("/ls7 {[600 159 200 150 200 150] 0 setdash} def\n")
ofp.write("/ls8 {[75 75] 0 setdash} def\n")
ofp.write("/ls9 {[600 150 75 150 75 150 75 150] 0 setdash} def\n")
ofp.write("/ls10 {[600 150 200 150 200 150 200 150] 0 setdash} def\n")
ofp.write("/L {lineto} def\n")
ofp.write("/m {moveto} def\n")
ofp.write("/g {setgray} def\n")
ofp.write("/s {stroke} def\n")

#define procedure to make filled rectangles
ofp.write("/F {\n")
ofp.write("   /blue exch def\n")
ofp.write("   /green exch def\n")
ofp.write("   /red exch def\n")
ofp.write("   /dy exch def\n")
ofp.write("   /dx exch def\n")
ofp.write("   /y exch def\n")
ofp.write("   /x exch def\n")
ofp.write("   s\n")
ofp.write("   red green blue setrgbcolor\n")
ofp.write("   x y m\n")
ofp.write("   0 dy rlineto\n")
ofp.write("   dx 0 rlineto\n")
ofp.write("   0 dy neg rlineto\n")
ofp.write("   closepath\n")
ofp.write("   fill\n")
ofp.write("} def\n")

ofp.write("%%EndProlog\n")
ofp.write("%%Page: 1 1\n")
ofp.write("gsave\n")

ofp.write("10 72 mul 32000 div 7.5 72 mul 24000 div scale\n")
if scale_option:
    ofp.write(f"{sfact} {sfact} scale\n")

#add the optional label for userid, dir, file, date and time
if id_option:
    ofp.write("gsave\n")
    ofp.write("32000 10 72 mul div 24000 7.5 72 mul div scale\n")
    ofp.write("0 2 m\n")
    ofp.write(f"({usname}          )show\n")
    ofp.write(f"({path}         )show\n")
    ofp.write("360 2 m (f001.sgf) stringwidth pop 2 div neg 0 rmoveto\n")
    ofp.write("(f001.sgf)show\n")
    ofp.write(f"720 2 m ({dattim}) stringwidth pop neg 0 rmoveto\n")
    ofp.write(f"({dattime})show\n")
    ofp.write("s\n")
    ofp.write("grestore\n")

ofp.write(f"{base_width} setlinewidth\n")

#execute this loop once for each input record
while not done:
    cbuf4=file_in.read(4)
    if not cbuf4:
        break
    buflen=int.from_bytes(bytes=cbuf4,byteorder="little")
    buflen*=2

    cbuf2=file_in.read(2*buflen)
    buffer=unpack(f"<{buflen}h",cbuf2)

    execute_buffer()

file_in.close()

#print the page
ofp.write("s\n")
ofp.write("showpage\n")
ofp.write("grestore\n")
ofp.write("%%Trailer\n")

ofp.close()
