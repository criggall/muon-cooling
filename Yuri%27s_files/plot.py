import numpy as np
from matplotlib import pyplot as plt

# RF period (325 MHz frequency):
T = 1/(325*10**6)*10**9 # ns

# Load data from output txt files:
initial_data = np.loadtxt('out1.txt')
final_data = np.loadtxt('out31.txt')

# Values for inital detector:
xi = []; yi = []; zi = []
pxi = []; pyi = []; pzi = []; ptotali = []
ti = []; mod_ti = []
for i in range(initial_data.shape[0]):
    xi.append(initial_data[i][0]/10) # mm -> cm
    yi.append(initial_data[i][1]/10)
    zi.append(initial_data[i][2]/10)
    px = initial_data[i][3]; py = initial_data[i][4]; pz = initial_data[i][5]
    pxi.append(px) # MeV/c
    pyi.append(py)
    pzi.append(pz)
    ptotali.append(np.sqrt(px**2+py**2+pz**2))
    t = initial_data[i][6]
    ti.append(t) # ns
    mod_ti.append(t % T)
    del px, py, pz, t

# Values for final detector:
xf = []; yf = []; zf = []
pxf = []; pyf = []; pzf = []; ptotalf = []
tf = []; mod_tf = []
for i in range(final_data.shape[0]):
    xf.append(final_data[i][0]/10) # mm -> cm
    yf.append(final_data[i][1]/10)
    zf.append(final_data[i][2]/10)
    px = final_data[i][3]; py = final_data[i][4]; pz = final_data[i][5]
    pxf.append(px) # MeV/c
    pyf.append(py)
    pzf.append(pz)
    ptotalf.append(np.sqrt(px**2+py**2+pz**2))
    t = final_data[i][6]
    tf.append(t) # ns
    mod_tf.append(t % T)
    del px, py, pz, t

# Plot px vs x:
plt.figure(1)
point_size = 5
plt.scatter(xi,pxi,color='blue',label='initial',s=point_size)
plt.scatter(xf,pxf,color='red',label='final',s=point_size)
plt.ylim(-100,100)
plt.xlabel('x (cm)')
plt.ylabel('p_x (MeV/c)')
plt.legend()
plt.title('p_x vs x')
plt.savefig('px_vs_x.png',dpi=300)

# Plot py vs y:
plt.figure(2)
plt.scatter(yi,pyi,color='blue',label='initial',s=point_size)
plt.scatter(yf,pyf,color='red',label='final',s=point_size)
plt.ylim(-100,100)
plt.xlabel('y (cm)')
plt.ylabel('p_y (MeV/c)')
plt.legend()
plt.title('p_y vs y')
plt.savefig('py_vs_y.png',dpi=300)

# Plot total p vs t:
plt.figure(3)
plt.scatter(mod_ti,ptotali,color='blue',label='initial',s=point_size)
plt.scatter(mod_tf,ptotalf,color='red',label='final',s=point_size)
# plt.ylim(0,400)
plt.xlabel('mod(t,T_RF)')
plt.ylabel('p_total (MeV/c)')
plt.legend()
plt.title('p_total vs mod(t,T_RF)')
plt.savefig('ptotal_vs_modtT.png',dpi=300)