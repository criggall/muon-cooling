import pandas as pd
import numpy as np
import math

# Read in ICOOL file:
icool_file = 'for003Pre6D.dat'
icool_df = pd.read_csv(icool_file)

# Reference particle data (not used in BLTrackFile):
ref = icool_df.iloc[0]

# Particle data (by event number):
data = icool_df.iloc[1:]

# Create txt file with column and unit labels:
f = open('output.txt','w')
f.write("#BLTrackFile\n")
f.write("#x y z Px Py Pz t PDGid EvNum TrkId Parent weight\n")
f.write("#cm cm cm MeV/c MeV/c MeV/c ns - - - - -\n")

for j in range(len(icool_df.values)-1):

    # Read ICOOL values by row:
    # print([float(x) for x in data.values[0][0].split()])
    valsICOOL = data.values[j][0].split()
    for k in range(len(valsICOOL)):
        valsICOOL[k] = float(valsICOOL[k])

    # Remove events with low event weight:
    eventWeight = valsICOOL[5]
    eventWeightOrder = math.floor(math.log(eventWeight, 10))
    valsBL = []
    if eventWeightOrder >= 0:

        # Rearrange and convert units to BLTrackFile format:
        valsBL.append(valsICOOL[6]*100) # x (m -> cm)
        valsBL.append(valsICOOL[7]*100) # y (m -> cm)
        valsBL.append(valsICOOL[8]*100) # z (m -> cm)
        valsBL.append(valsICOOL[9]*1000) # Px (GeV/c -> MeV/c)
        valsBL.append(valsICOOL[10]*1000) # Py (GeV/c -> MeV/c)
        valsBL.append(valsICOOL[11]*1000) # Pz (GeV/c -> MeV/c)
        valsBL.append(valsICOOL[4]*10**9) # t (s -> ns)
        valsBL.append(13.0) # PDGid (13 = muon)
        valsBL.append(valsICOOL[0]) # EvNum
        valsBL.append(-1.0) # TrkId (?)
        valsBL.append(-1.0) # Parent weight (try -1 for now)

    # Write to txt:
    valsStr = [str(x) for x in valsBL]
    f.write(' '.join(valsStr) + '\n')

f.close()

# Convert to .dat file:
with open('output.txt','r') as f:
    r = f.read()
with open('output.dat','w') as f:
    f.write(r)