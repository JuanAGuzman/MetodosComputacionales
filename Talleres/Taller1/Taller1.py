import os
import numpy as np
import matplotlib.pyplot as plt
import wget as w
from scipy.optimize import curve_fit

#Ejercicio 1

def barrido(f, i, n, p):
    signo = f(i)
    x = i
    ceros = []
    while x<=n:
        if (f(x)>0 and f(x+p)<0) or (f(x)<0 and f(x+p)>0):
            ceros.append(x)
            x+=p
        else:
            x+=p
    return ceros

m = 1
g = 9.8
y_0 = 0
v_0 = 50
C = 0.8
rho = 1.225
R = 0.05
A = np.pi * R**2

def gamma(C, rho, A, m):
    return np.sqrt(C*rho*A/2*m)

def y_i(t):
    return y_0 + v_0*t - g* t**2 /2

def v_i(t):
    return v_0 - g*t

def y_r(t):
    x = np.arctanh(gamma(C, rho, A, m)*v_0/np.sqrt(g))
    return y_0 + (np.log(np.cosh(x)/np.cosh(x - gamma(C, rho, A, m)*np.sqrt(g)*t)) / gamma(C, rho, A, m)**2)

def v_r(t):
    x = np.arctanh(gamma(C, rho, A, m)*v_0/np.sqrt(g))
    return np.sqrt(g)*np.tanh(x - gamma(C, rho, A, m)*np.sqrt(g)*t) / gamma(C, rho, A, m)

r_i = barrido(y_i, 1, 100, 0.001)
r_r = barrido(y_r, 1, 100, 0.001)

print('El tiempo de vuelo ideal es igual a ', r_i[-1])
print('El tiempo de vuelo resistivo es igual a ', r_r[-1])

x_1 = np.arange(0, r_i[-1], 0.001)
x_2 = np.arange(0, r_r[-1], 0.001)

y_01 = np.zeros(len(x_1))
y_02 = np.zeros(len(x_2))

y_1 = y_i(x_1)
v_1 = v_i(x_1)

y_2 = y_r(x_2)
v_2 = v_r(x_2)

plt.plot(x_1, y_1, label='Posicion (m)')
plt.plot(x_1, v_1, label='Velocidad (m/s)')
plt.xlabel('Tiempo (s)')
plt.legend()
plt.savefig('casoidealpregunta1.jpg')
plt.show()

def y2(t, C):
    x = np.arctanh(gamma(C, rho, A, m)*v_0/np.sqrt(g))
    return y_0 + (np.log(np.cosh(x)/np.cosh(x - gamma(C, rho, A, m)*np.sqrt(g)*t)) / gamma(C, rho, A, m)**2)

def barrido2(f, i, n, p, c):
    signo = f(i, c)
    x = i
    ceros = []
    while x<=n:
        if (f(x, c)>0 and f(x+p, c)<0) or (f(x, c)<0 and f(x+p, c)>0):
            ceros.append(x)
            x+=p
        else:
            x+=p
    return ceros

C = np.arange(0.1, 0.9, 0.1)
t = []
for i in C:
    t.append(barrido2(y2, 1, 100, 0.001, i))

plt.plot(C, t)
plt.xlabel('Coeficiente')
plt.ylabel('Tiempo de vuelo (s)')
plt.savefig('coeficiente_vuelo.jpg')
plt.show()

'''
#Ejercicio 2

def primo(x):
    for i in range(2, x):
        if x % i == 0:
            return False
    return True    

def primos(n):
    x = []
    for i in range(2, n, 1):
        if primo(i):
            x.append(i)
    return x

y = primos(1000)
x = np.arange(0, len(y), 1)

plt.scatter(x, y, s=1)
plt.xlabel('Posicion')
plt.ylabel('Numero')
plt.savefig('primos.jpg')
plt.show()


def f(n):
    a, b = 0, 1
    x = []
    while a<n:
        a, b = b, a+b
        x.append(a)
    return x

def aureo(n):
    return f(n)[-1]/f(n)[-2]

a_r = (1 + np.sqrt(5)) / 2


def e(x, r):
    return np.abs(x - r)/r

x = []
y = []
for i in range(2, 75):
    x.append(i)
    y.append(e(aureo(i), a_r))
    
plt.plot(x, y, c='r')
plt.scatter(x, y, c='b', s=10)
plt.xlabel('Numeros generados')
plt.ylabel('Error')
plt.savefig('error_aureo.jpg')
plt.show()


#Ejercicio 3

url = 'https://raw.githubusercontent.com/asegura4488/MetodosCompu2021/main/Week2/Data/Hw_data.dat'
w.download(url, str(os.getcwd())+'/Hw_data.dat')
archivo = open(str(os.getcwd())+'/Hw_data.dat', 'r')
f = archivo.readlines()

x, y = [], []
for i in range(len(f)):
    x.append(float(f[i].split(' ')[0]))
    y.append(float(f[i].split(' ')[1]))

plt.plot(x, y)
plt.savefig('ejercicio3.jpg')

def f(t, A, B):
    return A/ (1 + np.exp(-B*t) )

popt, pcov = curve_fit(f, x, y)

A, B = popt[0], popt[1]

print('El valor de A es igual a ',A, ' y el de B ', B)


#Ejercicio 4

def f(x):
    return (3* x**5) + (5* x**4) - (x**3)

def newton(f, x, n):
    h = 0.00000001
    a = n
    y = f(x)
    while a>0:
        if y==0:
            return x
        a-=1
        x = x - f(x)/(f(x+h)-f(x-h)/(2*h))
        y = f(x)
    return x

def barrido(f, i, n, p):
    signo = f(i)
    x = i
    ceros = []
    while x<=n:
        if (f(x)>0 and f(x+p)<0) or (f(x)<0 and f(x+p)>0):
            ceros.append(x)
            x+=p
        else:
            x+=p
    return [newton(f, x, n) for x in ceros]

print('Las raices se encuentran en ', barrido(f, -100, 100, 0.001))

#Ejercicio 6

def f(x):
   return np.exp(-x**2)
def simpsons38(f,a,b):
    m1=(2*a+b)/3
    m2=(a+2*b)/3
    integral=(b-a)/8*(f(a)+3*f(m1)+3*f(m2)+f(b))
    return integral
def simpsons13(f,a,b):
    m=(a+b)/2
    integral=(b-a)/6*(f(a)+4*f(m)+f(b))
    return integral
prin(simpsons38)
a=0
b=1
n=100
h=(b-a)/n
suma=0
for i in range(n):
    b=a+h
    area=simpsons38(f,a,b)
    suma=suma+area
    a=b
print(suma)
vt=0.746824132812427
error1=abs((vt-suma)/ vt)*100
print("Error% =", error1)
a=0
b=1
n=100
hr=(b-a)/n
sumat=0
for i in range(n):
    b=a+h
    areat=simpsons13(f,a,b)
    sumat=sumat +areat
    a=b
print(sumat)
errort=abs((vt-sumat)/ vt)*100
print("Error% =", errort)
'''