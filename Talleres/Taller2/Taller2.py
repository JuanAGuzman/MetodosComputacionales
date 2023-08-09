import os
import numpy as np
import wget as w
import matplotlib.pyplot as plt
from matplotlib import rc
from tqdm import tqdm
from mpl_toolkits.mplot3d import axes3d


#Pregunta 1

class MyRandom():
    def __init__(self, seed = 15, method='simple'):
        self.r = seed
        self.method = method
        if method=='simple':
            self.a = 57
            self.c = 1
            self.M = 265
        elif method == 'drand48':
            self.a = int('5DEECE66D',16)
            self.c = int('B',16)
            self.M = 2**48
        else:
            print('Generador no reconocido')
    def Random(self):
        r = (self.a*self.r + self.c)%self.M
        self.r = r
        return r/float(self.M)
    def TestMethod(self, Npoints, moment, seed_ = 32, method_ = 'simple'):
        rand = MyRandom(seed = seed_, method = method_)
        array = np.zeros(Npoints)
        for i in range(Npoints):
            array[i] = rand.Random()
        return np.sqrt(Npoints)* np.abs(  np.mean(array**moment) - 1./(1.+moment) )

def FillPoints(seed_, method_, Npoints):
    rand = MyRandom(seed = seed_, method = method_)
    points = np.zeros(Npoints)
    for i in tqdm(range(Npoints)):
        points[i] = rand.Random()
    return points

Npoints = 10**3
Nsimple = FillPoints(165, 'simple', Npoints)
Nrand48 = FillPoints(695, 'drand48', Npoints)

def evaluador(x, k):
    return sum([x[i]*x[i+int(k)] for i in range(len(x)-int(k))])/len(x)

k = 20
print(evaluador(Nsimple, k))
print(evaluador(Nrand48, k))

x = np.arange(1, 100, 1)
y = [evaluador(Nsimple, i) for i in x]
plt.plot(x, y)
plt.show()

x = np.arange(1, 100, 1)
y = [evaluador(Nrand48, i) for i in x]
plt.plot(x, y)
plt.show()


#Pregunta 2

def vol(d, r, n):
    x = np.random.uniform(0, 1, tuple([d, n]))*2 - 1
    y = r*np.array([i[np.where(sum(x**2)<1)] for i in x])
    return r*x, y

cubo, esfera = vol(3, 2, 2000)
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1,projection='3d')
ax.view_init(10,60)
ax.scatter(esfera[0], esfera[1], esfera[2])
plt.show()

#Pregunta 3

cubo, esfera = vol(3, 1, 200)
def f(x):
    y = [np.exp(np.sqrt(np.sum(np.power(x[:,i], 2)))) for i in range(len(x[0]))]
    return y

def Integral(esfera, Npoints):
    Integral = 4*np.pi*np.average([f(esfera) for i in range(Npoints)])/3
    return Integral


x = np.arange(10, 5010, 100)
x_0 = np.array([Integral(vol(3, 1, i)[1], 100) for i in x])
y = x_0/(4*np.pi*(np.e-2))

plt.plot(x, y)
plt.show()


#Pregunta 4

def f(x):
    r = np.sqrt(np.sum(np.power(x, 2)))
    return np.cos(r)*np.exp(-np.power(r, 2))

def g(x):
    return x[0]

def integral_metropolis(f,g,dimension, steps):
    v_gx = np.zeros(steps)
    old_x = np.ones(dimension)
    old_prob = f(old_x)
    for i in range(steps):
        new_x = np.random.normal(old_x, 0.5,dimension)
        new_prob = f(new_x)
        aceptacion = new_prob / old_prob
        if aceptacion > np.random.rand():
            old_x = new_x
            old_prob = new_prob
        v_gx[i] = g(old_x)
    return 2*np.std(v_gx)/np.sqrt(steps)

print(integral_metropolis(f, g, 2, 10**5), integral_metropolis(f, g, 3, 10**5))

#Pregunta 5

url = 'https://raw.githubusercontent.com/asegura4488/MetodosCompu2021/main/Week7/data/MCMC_data.dat'
w.download(url, str(os.getcwd())+'/pregunta5_data.dat')
archivo = open(str(os.getcwd())+'/pregunta5_data.dat', 'r')
f = archivo.readlines()

x, y = [], []
for i in range(len(f)):
    x.append(float(f[i].split(' ')[0]))
    y.append(float(f[i].split(' ')[1]))
x = np.array(x)
y = np.array(y)

def f_1(x, c):
    return c[0]*x + c[1]

def f_2(x, c):
    return c[0]* x**2 + c[1]*x + c[2]

def f_3(x, c):
    return c[0]*np.exp(c[1]*x)

def metropolis_fit(t,D,model,dim,n):
    muestras = np.zeros([dim,n])
    old_theta = np.ones(dim)
    old_prob = np.exp(-np.sum((D - model(t,old_theta))**2)/2)
    for i in range(n):
        new_theta = old_theta + np.random.normal(0, 0.5, dim)
        new_prob = np.exp(-np.sum((D - model(t,new_theta))**2)/2)
        aceptacion = new_prob / old_prob
        if aceptacion > np.random.rand():
            old_theta = new_theta
            old_prob = new_prob
        muestras[:,i] = old_theta
    return muestras

def calcular_parametro(samples,inicio,dimension):
    return np.mean(samples[dimension,inicio:])* 10**5

def mostrar_ajuste_1(t,D,samples,inicio):
    c0 = calcular_parametro(samples,inicio,0)
    c1 = calcular_parametro(samples,inicio,1)
    print(c0,c1)
    y = c0*t + c1
    plt.plot(t,D,'o')
    plt.plot(t,y,c='r')
    plt.show()

def mostrar_ajuste_2(t,D,samples,inicio):
    c0 = calcular_parametro(samples,inicio,0)
    c1 = calcular_parametro(samples,inicio,1)
    c2 = calcular_parametro(samples,inicio,2)
    print(c0,c1,c2)
    y = c0* t**2 + c1*t + c2
    plt.plot(t,D,'o')
    plt.plot(t,y,c='r')
    plt.show()

def mostrar_ajuste_3(t,D,samples,inicio):
    c0 = calcular_parametro(samples,inicio,0)
    c1 = calcular_parametro(samples,inicio,1)
    print(c0,c1)
    y = c0*np.exp(c1*t)
    plt.plot(t,D,'o')
    plt.plot(t,y,c='r')
    plt.show()

samples_1 = metropolis_fit(x, y, f_1, 2, 100000)
samples_2 = metropolis_fit(x, y, f_2, 3, 100000)
samples_3 = metropolis_fit(x, y, f_3, 2, 100000)
mostrar_ajuste_1(x, y, samples_1, 1)
mostrar_ajuste_2(x, y, samples_2, 1)
mostrar_ajuste_3(x, y, samples_3, 1)

#Pregunta 6

url = 'https://raw.githubusercontent.com/asegura4488/MetodosCompu2021/main/Week7/data/Likelihood.dat'
w.download(url, str(os.getcwd())+'/pregunta6_data.dat')
archivo = open(str(os.getcwd())+'/pregunta6_data.dat', 'r')
f = archivo.readlines()

x = np.array([float(i) for i in f])

def likelihood(x,mu,sigma):
    return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def Chi2(f, x, mu, sigma):
    x_fit = np.log(f(x, mu, sigma))
    chi2 = []
    for i in range(len(x)):
        chi2.append((x[i]-x_fit[i])**2/(2*sigma**2))
    return np.array(chi2)

mu = np.linspace(0, 10, 10)
sigma = np.linspace(0, 10, 10)
z = Chi2(likelihood, x, mu, sigma)

fig, ax = plt.subplots()
c = ax.contour(mu, sigma, z, 20)
ax.clabel(c, inline=1, fontsize=10)
plt.show()

ind = np.where(z == np.min(z))

zmin = z[ind]

Mumin, Sigmamin = mu[ind[0]], sigma[ind[1]]

print("El valor aproximado de Mu es:", Mumin,"y el valor aproximado de Sigma es:",Sigmamin)

