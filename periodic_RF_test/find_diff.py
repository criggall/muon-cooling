import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import sys
import os
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)

from functions.read_g4bl_data import readDetData
from functions.set_plot_settings import setPlotSettings

setPlotSettings(font=True)

########################################################################

data = readDetData('zntuple.txt')
# data = readDetData('zntuple2.txt')
# data = readDetData('zntuple3.txt')

# data = data[data['EventID']==1]

# remove last period (for fringe-field mitigation)
# 2 periods HFOFO + 1 at start + 1 at end
period_len = 4200
data = data[data['z'] <= 500+period_len*2]

data0 = data[data['z']==500] # start
data1 = data[data['z']==500+period_len*2] # after 2 turns (periods)

x0 = data0['x'].values
y0 = data0['y'].values
x1 = data1['x'].values
y1 = data1['y'].values

diffs = []
for i in range(len(x0)):
    diffs.append(np.sqrt((x1[i]-x0[i])**2+(y1[i]-y0[i])**2))
    if diffs[i] <= np.min(diffs):
        min_index = i

print(min(diffs))

########################################################################

# min_soln = {
#     'x' : data0['x'].values[min_index],
#     'y' : data0['y'].values[min_index],
#     'z' : data0['z'].values[min_index],
#     'xp' : data0['px'].values[min_index]/data0['pz'].values[min_index],
#     'yp' : data0['py'].values[min_index]/data0['pz'].values[min_index],
#     't' : data0['t'].values[min_index]
# }

# print(min_soln)

# diff = np.sqrt((x1-x0)**2+(y1-y0)**2)
# print(diff)

# plt.plot(data['x'],data['y'])
# plt.scatter(x0,y0, color='tab:green', label='intial')
# plt.scatter(x1,y1, color='tab:red', label='final')
# plt.xlabel('$x$ [mm]')
# plt.ylabel('$y$ [mm]')
# plt.legend(loc='lower right', fontsize=10)
# plt.title(f'Difference = {diffs} mm')
# plt.savefig('diff.png', dpi=300)