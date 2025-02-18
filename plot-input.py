# Import libraries:
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

dir = 'Yuri-David-File/'
file = dir+'initial.dat'

df = pd.read_csv(file, comment='#', delim_whitespace=True, 
                 names=["x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackID", "ParentID", "Weight"])

px = df["Px"].to_numpy()
py = df["Py"].to_numpy()
pz = df["Pz"].to_numpy()

ptotal = []
for i in range(len(px)):
    ptotal.append(np.sqrt(px[i]**2+py[i]**2+pz[i]**2))

plt.hist(ptotal,color='green',range=(0,400),bins=80) # bins of 5
# plt.hist(ptotal,color='green',bins=80) # bins of 5
plt.xlabel('p_total (MeV/c)',loc='right')
plt.ylabel('Events',loc='top')
plt.savefig(dir+'input_ptotal_distribution.png',dpi=300)