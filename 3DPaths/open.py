import matplotlib.pyplot as plt
import numpy as np

x = [line.strip('\n').split(',') for line in open("data/x.csv")]
x = np.array([(float(a),float(b)) for a,b in x])
y = [line.strip('\n').split(',') for line in open("data/y.csv")]
y = np.array([(float(a),float(b)) for a,b in y])

#fig,ax = plt.subplots()

plt.scatter(x[:,0],x[:,1],color="red",marker="o",label="x",s=20)
plt.scatter(y[:,0],y[:,1],color="blue",marker="o",label="y",s=20)

z = np.linspace(1,400,100)
coefx = np.polyfit(x[:,0],x[:,1],15)
poly1d_x = np.poly1d(coefx) 
plt.plot(z,poly1d_x(z), color="red")
coefy = np.polyfit(y[:,0],y[:,1],15)
poly1d_y = np.poly1d(coefy) 
plt.plot(z,poly1d_y(z), color="blue")

plt.clf()

ax = plt.axes(projection ='3d')
ax.plot3D(z, poly1d_x(z), poly1d_y(z), 'green')
ax.set_title('3D Trajectory')

plt.clf() 
period = 10
x = np.concatenate([poly1d_x(z) for i in range(period)])
y = np.concatenate([poly1d_y(z) for i in range(period)])

z = np.linspace(0,400*period,100*period)
ax = plt.axes(projection ='3d')
ax.plot3D(z, x/z**0.2, y/z**0.2, 'green')
ax.plot3D(z, 0.7*x/z**0.2, 0.7*y/z**0.2, 'blue')
ax.plot3D(z, 0.4*x/z**0.2, 0.4*y/z**0.2, 'red')
ax.set_title('3D Trajectory')

plt.show()
