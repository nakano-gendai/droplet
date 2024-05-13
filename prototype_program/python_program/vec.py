import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi, sqrt

ext = 'pdf' # 保存ファイルの拡張子　pdf,svg,pngなど
datadir =""
dfile1   = datadir + '0_1000.d'

# x = np.loadtxt(dfile1, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# y = np.loadtxt(dfile1, usecols = 1, dtype = 'float64')
# vx = np.loadtxt(dfile1, usecols = 2, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# vy = np.loadtxt(dfile1, usecols = 3, dtype = 'float64')

# fig = plt.figure(figsize = (2.4, 2.4), dpi = 100, linewidth = 0)
# ax = fig.add_subplot(111)
# ax.grid()

# ax.quiver(x, y, vx, vy, color = "black", width=0.001,
#           angles = 'xy', scale_units = 'xy', scale = 1.4)
# plt.streamplot(x, y, vx, vy)

# plt.show()


n = 1000
x = np.linspace(0, 2*pi, n)
y = np.linspace(0, 2*pi, n)
X, Y = np.meshgrid(x, y)

# Taylor-Green Vortices
u = sin(X)*cos(Y)
v = - cos(X)*sin(Y)
speed = sqrt(u**2 + v**2)

# Plot
plt.figure(1)
plt.streamplot(X, Y, u, v, density=3, color='k', arrowstyle='-', linewidth=0.6)
plt.contourf(X, Y, speed, 100, cmap='viridis')
plt.show()