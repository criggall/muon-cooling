from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
 
fig = plt.figure()


z = np.linspace(1, 5*np.pi, 1000)
phase1 = 0
phase2 = 0
x = np.cos(5 * z + phase1)/z**0.5 #+ z**3*0.0001*np.sin(7001001 * z + phase2)  +  np.sin(2 * z + phase2) 
y = np.sin(5 * z + phase1)/z**0.5 #+ z**3*0.0001*np.cos(7001001 * z + phase2)  +  np.cos(2 * z + phase2)
 
# plotting
ax = plt.axes(projection ='3d')
ax.plot3D(z, x, y, 'green')
ax.plot3D(z, 0.4*x, 0.4*y, 'blue')
ax.set_title('3D Trajectory')

# ax = plt.axes()
# ax.plot(z, x, color="red")
# ax.plot(z, y, color="blue")
plt.show()
