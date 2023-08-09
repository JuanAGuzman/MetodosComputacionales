import numpy as np
import matplotlib.pyplot as plt

Data=np.array([10,20,5])
bkg=np.array([6,15,10])
Signal=np.array([2,10,3])

def likelihood(s, b, N, c):
    return np.prod([(np.exp(-(c[0]*s + c[1]*b)) * (c[0]*s + c[1]*b)**i *np.exp(-(c[1]-1)**2 / (2* 0.1**2))) / (np.math.factorial(i)*np.sqrt(2*np.pi* 0.1**2)) for i in range(N)])

def metropolis_fit(s, b, N, f, dim, n):
    muestras = np.zeros([dim,n])
    old_theta = np.ones(dim)
    old_prob = f(s, b, N, old_theta)
    for i in range(n):
        new_theta = old_theta + np.random.normal(0, 0.5, dim)
        new_prob = f(s, b, N, new_theta)
        aceptacion = new_prob / old_prob
        if aceptacion > np.random.rand():
            old_theta = np.abs(new_theta)
            old_prob = new_prob
        muestras[:,i] = old_theta
    return muestras

def calcular_parametros(samples,inicio,dimensiones):
    return [np.mean(samples[i,inicio:]) for i in range(dimensiones)]

muestras = [metropolis_fit(Signal[i], bkg[i], Data[i],likelihood, 2, 10000) for i in range(len(Data))]

x = np.transpose([calcular_parametros(muestras[0], 10, 2) for i in range(len(Data))])
mu = np.mean(x[0])
epsilon = np.mean(x[1])
print('El valor de mu que se obtuvo fue ', mu, ' y el de epsilon ', epsilon)


def f(x):
    return 2**(-7) *np.sum(np.power(x, 1), axis=1)**2

def Integral(f, a, b, dim, n):
    x = np.random.uniform(a, b, (n, dim))
    return np.sum(f(x))* (b-a)**n/ n

print('El valor de la integral que se obtuvo fue ', Integral(f, 0, 1, 8, 1000000))    