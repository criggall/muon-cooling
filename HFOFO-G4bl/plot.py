import numpy as np
from matplotlib import pyplot as plt

# Load data from output txt files:
initial_data = np.loadtxt('out3.txt') # out1.txt and out2.txt empty - use out3.txt for now
final_data = np.loadtxt('out31.txt')

# Values for inital detector:
xi = []; yi = []; zi = []
pxi = []; pyi = []; pzi = []; ptotali = []
ti = []
for i in range(initial_data.shape[0]):
    xi.append(initial_data[i][0]/10) # mm -> cm
    yi.append(initial_data[i][1]/10)
    zi.append(initial_data[i][2]/10)
    pxi.append(initial_data[i][3]) # MeV/c
    pyi.append(initial_data[i][4])
    pzi.append(initial_data[i][5])
    ptotali.append(initial_data[i][3]+initial_data[i][4]+initial_data[i][5])
    ti.append(initial_data[i][6]) # ns

# Values for final detector:
xf = []; yf = []; zf = []
pxf = []; pyf = []; pzf = []; ptotalf = []
tf = []
for i in range(final_data.shape[0]):
    xf.append(final_data[i][0]/10) # mm -> cm
    yf.append(final_data[i][1]/10)
    zf.append(final_data[i][2]/10)
    pxf.append(final_data[i][3]) # MeV/c
    pyf.append(final_data[i][4])
    pzf.append(final_data[i][5])
    ptotalf.append(final_data[i][3]+final_data[i][4]+final_data[i][5])
    tf.append(final_data[i][6]) # ns

# Plot px vs x:
plt.figure(1)
plt.scatter(pxi,xi,color='blue',label='initial')
plt.scatter(pxf,xf,color='red',label='final')
plt.xlabel('x (cm)')
plt.ylabel('p_x (MeV/c)')
plt.legend()
plt.title('p_x vs x')
plt.savefig('px_vs_x.png',dpi=300)

# Plot py vs y:
plt.figure(2)
plt.scatter(pyi,yi,color='blue',label='initial')
plt.scatter(pyf,yf,color='red',label='final')
plt.xlabel('y (cm)')
plt.ylabel('p_y (MeV/c)')
plt.legend()
plt.title('p_y vs y')
plt.savefig('py_vs_y.png',dpi=300)

# Plot total p vs t:
plt.figure(3)
plt.scatter(ptotali,ti,color='blue',label='initial')
plt.scatter(ptotalf,tf,color='red',label='final')
plt.xlabel('t (ns)')
plt.ylabel('p_total (MeV/c)')
plt.legend()
plt.title('p_total vs t')
plt.savefig('ptotal_vs_t.png',dpi=300)

#######################################################

# x_vals = []; y_vals = []; z_vals = []
# px_vals = []; py_vals = []; pz_vals = []
# t_vals = []
# for j in range(31):

#     # Load data from output txt files:
#     data = np.loadtxt('out'+str(j+1)+'.txt')

#     # Values for each detector:
#     x = []; y = []; z = []
#     px = []; py = []; pz = []
#     t = []
#     for i in range(data.shape[0]):
#         x.append(data[i][0]/10) # mm -> cm
#         y.append(data[i][1]/10)
#         z.append(data[i][2]/10)
#         px.append(data[i][3]) # MeV/c
#         py.append(data[i][4])
#         pz.append(data[i][5])
#         t.append(data[i][6]) # ns

#     # Values for entire channel:
#     x_vals.append(x)
#     y_vals.append(y)
#     z_vals.append(z)
#     px_vals.append(px)
#     py_vals.append(py)
#     pz_vals.append(pz)
#     t_vals.append(t)