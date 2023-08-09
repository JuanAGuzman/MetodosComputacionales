import os
import wget as w
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import *

# Colision de dos discos rigidos

m = 1.5
u_1= 2
e = 0.8
k = 0.5
theta = np.arange(0, np.pi/2, np.pi/2)

def solver(x):
    return np.matmul(np.linalg.inv(np.arange([[np.sin(x), np.cos(x), -np.sin(x), -np.cos(x), 1, 1],
                                              [0, m, 0, 1, 0, 0],
                                              [-np.cos(x), np.sin(x), np.cos(x), -np.sin(x), 0, 0],
                                              [np.sin(x), np.cos(x), 0, 0, -k, 0],
                                              [0, 0, np.sin(x), np.cos(x), 0, k]])), np.arange(0, 0, m*u_1, e*u_1*np.cos(x), u_1*np.sin(x), 0))

sol = [solver(i) for i in theta]
print(sol)


# Precesion de la orbita de Mercurio


dt = 0.01
alpha = 1.1 * 10**(-8)
e = 0.205630
a = 0.387098
T = 88
c = a*(1+e)
G = 2.9597 * 10**(-4)
t = np.arange(0, T, dt)
r = np.zeros((2, len(t)))
v = np.zeros((2, len(t)))

r[0][0], r[1][0] = c, 0
v[0][0], v[1][0] = 0, np.sqrt(G*(1-e)/c)

def a(r_x, r_y):
    r = r_x**2 + r_y**2
    return -(G/ r)*(1 + (alpha/ r))

for i in range(len(t)-1):
    r[0][i+1] = r[0][i] + v[0][i]*dt + (a(r[0][i], r[1][i])*(dt)**2 * np.cos(np.tan(r[1][i]/r[0][i]))/ 2)
    v[0][i+1] = v[0][i] + (a(r[0][i+1], r[1][i+1])+a(r[0][i], r[1][i]))*dt*np.cos(np.tan(r[1][i]/r[0][i]))/2
    r[1][i+1] = r[1][i] + v[1][i]*dt + (a(r[0][i], r[1][i])*(dt)**2 * np.sin(np.tan(r[1][i]/r[0][i]))/ 2)
    v[1][i+1] = v[1][i] + (a(r[0][i+1], r[1][i+1])+a(r[0][i], r[1][i]))*dt*np.sin(np.tan(r[1][i]/r[0][i]))/2

x = [i[0] for i in r]
y = [i[1] for i in r]


# Tensor de Inercia de una estrella naciente

url = 'https://raw.githubusercontent.com/asegura4488/MetodosCompu2021/main/Week10/data/CuerposCelestes.dat'
w.download(url, str(os.getcwd())+'/estrella_data.dat')
archivo = open(str(os.getcwd())+'/estrella_data.dat', 'r')
f = archivo.readlines()

E = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
r = np.array([[float(i.split(' ')[0]), float(i.split(' ')[1]), float(i.split(' ')[2])] for i in f])
I = sum([(np.dot(i, i)*E - np.kron(i, i).reshape((3, 3))) for i in r])
print(I)

valores, vectores = np.linalg.eig(I)
print(valores)
print(vectores)

x = [i[0] for i in r]
y = [i[1] for i in r]
z = [i[2] for i in r]

fig = plt.figure()
ax = Axes3D(fig, elev = 18, azim = 9)
ax.scatter(x, y, z, marker='.')
plt.show()

# Series de Fourier



# Transformada de Fourier

url = 'https://raw.githubusercontent.com/asegura4488/MetodosCompu2021/main/Week10/data/ManchasSolares.dat'
w.download(url, str(os.getcwd())+'/manchas_data.dat')
archivo = open(str(os.getcwd())+'/manchas_data.dat', 'r')
f = archivo.readlines()

print(f[0].split('  ')[0])
data = np.array([float(j) for i in f for j in i.split() if j.isdigit()])
#data = np.array([i for i in data if i[0]>=1900])
print(data)

# Redes Neuronales


