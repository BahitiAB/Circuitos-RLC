import numpy as np
import matplotlib.pyplot as plt



##################################################
##  Newton rapshon
##################################################
##################################################
import numpy as np 
from scipy.misc import derivative
import matplotlib.pyplot as plt


def newton_raphson(f, x0, ERROR, MAX_ITER):
    """
        Calcula un solucion de la ecuacion no lineal f(x)=0
        input
            f: funcion 
            x0: punto de inicio

        ouput
            x: solucion de la ecuacion no lineal f(x)=0
    """
#    xk = np.array([x0])
    for i in range(MAX_ITER):
        x = x0 - f(x0)/derivative(f, x0)
        #xk = np.append(xk, x)
        e = np.abs(x-x0)
        print("Iter{}: x = {}\tERROR = {}".format(i, x, e))
        if e < ERROR:
            return x
        x0 = x
    return x





## resonancias

def resonancias(a, b, f, ERROR=1e-5):
    x0 = (a+b)/2
    sol = newton_raphson(f, x0, ERROR, MAX_ITER=100)
    print("Solucion de la ecuacion f(x)=0\nx = {}\tf({}) = {}".format(sol, sol, f(sol)))
    # Grafica de la funcion x0
    DISC = 1000
    x = np.linspace(a, b, DISC)
    plt.plot(x, f(x), label="f(x)")
    plt.plot(sol, f(sol), "o", label="solucion")
    plt.grid()
    plt.legend()
    plt.show()


    

def inicializacion():

    print("Ingrese intervalo [a, b] a analizar")
    a0 = float(input("Ingrese a: "))
    b0 = float(input("Ingrese b: "))
    
    R1 = 1e-2
    R2 = 1e-2
    R3 = 1e-2

    L1 = 18.31e-3    # Asignar valor real del inductancia
    L2 = 11.27e-3    # Asignar valor real del inductancia
    L3 = 17.51e-3    # Asignar valor real del inductancia

    C1 = 1e-2    # Asignar valor real de la capacitancia
    C2 = 1e-2    # Asignar valor real de la capacitancia
    C3 = 1e-2    # Asignar valor real de la capacitancia



    a = lambda w: R3**2 +(L3*w - 1/(C2*w)- 1/(C3*w))**2
    b = lambda w: (R2 + R3/(a(w)*w**2*C2**2))**2 + (L2*w + L3/(a(w)*w*C2**2) + 2*L3/(a(w)*w*C2*C3) - w*L3**2/(a(w)*C2) - R3**2/(a(w)*w*C2)-1/(w*C1)-1/(a(w)*w**3*C2**2*C3)-1/(a(w)*w**3*C2*C3**2))**2


    m  = (-1)*(C2+C3)**2/(C1*C2**4*C3**4)
    n = -(C2+C3)/(C1**2*C2**2*C3**2) + (6*C2*C3+4*C2**2+2*C3**2)*L3/(C1*C2**4*C3**2) - 2*(C2+C3)*R3**2 - R3**2/(C1*C2**4)
    l = lambda w: a(w)*(2*(C2+C3)*L2/(C1*C2**2*C3**2) + (2*C2+C3)*L3/C1**2*C2**2*C3 - R3**2/(C1**2*C2) - 2*R2*R3/(C1*C2**2)) +\
            -(6*C2*C3+6*C2**2+C3**2)*L3**2/(C1*C2**4*C3**2) - R3**4/(C1*C2**2) + 2*(C2+C3)*L3*R3**2/(C1*C2**3*C3**2)
    p = lambda w: a(w)**2*(L2/C1**2 - R2**2/C1) + \
            ((2*L2*R3**2)/(C1*C2) -L3**2/(C1**2*C2) - 2*(2*C2+C3)*L2*L3/(C1*C2**2*C3))*a(w) +\
            2*(2*C2+C3)*L3**3/(C1*C2**2*C3) - 2*L3**2*R3**2/(C1*C2**2)

    r = lambda w: L1*a(w)**2*b(w) - L2**2/C1*a(w)**2 + 2*L2*L3**2/(C1*C2)*a(w) - L3**4/(C1*C2**2)

    f = lambda w: m + w**2*n + l(w)*w**4 + p(w)*w**6 + r(w)*w**8

    return a0, b0,  np.vectorize(f)

    
##################################################



 
"""
 
w = np.linspace(1, 50, 100)
fvec = np.vectorize(f)
 
plt.plot(w, fvec(w))
plt.grid()
plt.show()


#fvec = np.vectorize(f)
a0, b0, f = inicializacion()
resonancias(a0, b0, f)
"""

a0, b0 , f = inicializacion()
x = np.linspace(a0, b0, 100)
plt.plot(x, f(x))
plt.show()
