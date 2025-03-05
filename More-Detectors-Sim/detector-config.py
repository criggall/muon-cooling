# Channel length:
channel = 130200 # 31 periods of length 4200

# Number of detectors:
num = 1302

# Detector spacing:
spacing = channel/num

# Create txt file for detector configurations:
f = open('detector-config.txt','w')

for i in range(num):
    pos = spacing*i
    f.write("place Det"+str(i+1)+" z=-9.8+$zshift+"+str(pos)+"\n")

# Create txt file for detector output:
f = open('detectors.txt','w')

for i in range(num):
    f.write("virtualdetector Det"+str(i+1)+" file=out"+str(i+1)+" format=ascii radius=360 color=0,1,0 length=0.001 material=Vacuum\n")