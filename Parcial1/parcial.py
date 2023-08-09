import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from numpy.polynomial.legendre import leggauss
import sympy as sym

def derivada1(f, x, h):
    return (f(x+h) - f(x-h))/ (2*h)

def derivada2(f, x, h):
    return (f(x+h) - 2*f(x) + f(x-h))/ h**2

def ts(f, a, b):
    return (b-a)*(f(a) + f(b))/2

def tc(f, a, b, n):
    h = (b-a)/n
    return (h/2)*(f(a)+f(b)) + h*sum(f(np.linspace(a, b, n)))

def ss(f, a, b):
    h = (b-a)/2
    return h*(f(a) + 4*f(h) + f(b))/3

def sc(f, a, b, n):
    h =(b-a)/n
    return (h/3) * (f(a) + 4*sum(f(np.linspace(a, b, n))[0::2]) + 2*sum(f(np.linspace(a, b, n))[1::2]) + f(b))


xi, xf, Npoints = -1.,1.,1000
h = (xf-xi)/float(Npoints)

def NewtonMethod(f,df,xn,error,it,precision=0.001,iterations=1000):
    
    h_ = 1.0e-4
    
    while error > precision and it < iterations:
        
        try:
            
            xn1 = xn - f(xn)/df(f,xn,h_)+1.0e-10
            error = np.abs((xn1-xn)/xn1)        
            
        except ZeroDivisionError:
            print("Division by Zero")
            
        xn = xn1
        it += 1
    
    #print(it)
    if it == iterations:
        return False
    else:
        return xn1

root = NewtonMethod(Function,Derivative,1,100,0)
def GetRoots(f,df, X, precision_=0.001, tolerancia=5):
    
    Roots = []
    
    for i in X:
        
        root = NewtonMethod(f,df,i,100,0,precision=precision_)
        
        if root != False:
            if round(root,tolerancia) not in Roots:
                Roots.append(round(root,tolerancia))
            
      
    return Roots

def CreateLegendPoly(n):
    x1 = sym.Symbol('x', real=True)
    y = sym.Symbol('y', real=True)
    
    y = (x1**2-1)**n
    poly = sym.diff(y,x1,n)/( 2**n * math.factorial(n) )
    
    return poly


def GetWeight(f,df,xk):
    return 2./( (1-xk**2)*(df(xk))**2 )

def Getleggauss(n):
    
    if n == 0 or n == 1:
        return 0,0
    
    Legendre = []
    DLegendre = []
    Weights = []
    
    x = sym.Symbol('x', real=True)  
    
    for i in range(0,n+1):
        
        poly = CreateLegendPoly(i)
        Legendre.append(poly)
        DLegendre.append( sym.diff(poly,x,1) )
    
     
    xi = np.linspace(-1,1,200)
    
    pn = sym.lambdify([x], Legendre[n],'numpy')
    dpn = sym.lambdify([x], DLegendre[n],'numpy')

    Roots = GetRoots(pn,Derivative, xi, 0.00001,tolerancia=8)
    Roots.sort()
    
    
    for j in Roots:
        Weights.append(round(GetWeight(pn,dpn,j),8))
        
    Roots = np.array(Roots)
    Weights = np.array(Weights)
        
    return Roots, Weights

deg = 10
Roots, Weights = Getleggauss(deg)

a = 0
b = 0.25*np.pi

f = lambda x : np.sin(x)

t = 0.5*( (b-a)*Roots + a + b  )
Integral = 0.5*(b-a)*sum( Weights*f(t) )