import numpy as np
import matplotlib.pyplot as plt
import sys


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
    if(abs(f(x0)) < ERROR):
            return x0
    
    #DELTA = 1e-3
    for i in range(MAX_ITER):   
        x = x0 - f(x0)/derivative(f, x0)
        #xk = np.append(xk, x)
        e = abs(f(x))
        #step = abs(x-x0)
#        print("Iter{}: x = {}\tERROR = {}".format(i, x, e))
        #if step < DELTA or e < ERROR:
        if e <ERROR:
            return x
        x0 = x
    return x




def biseccion(f, a, b, ERROR, NUM_ITERATIONS):
    """
        Calcula algún cero de la función f en el intervalo [a, b].
    """
#    DELTA = 1e-8
    if f(a)*f(b) < 0:
        x = (a+b)/2
        k = 0
        while  abs(f(x)) > ERROR and k < NUM_ITERATIONS:
            if f(x)*f(a) < 0:
                b = x
            else:
                a = x

            x = (a+b)/2
            k +=1

            print("Iter {}: a={}\tb={}\tx= {}\tf(x) = {}".format(k, a, b, x, f(x)))
            
        if k == NUM_ITERATIONS:
            print("Maximo numero de iteraciones alcanzadas(biseccion)")
            return x
        else:
            print("Maxima cota error superada(biseccion)")
            return x

    else:
        print("Aplique otro método más adecuado.")




## resonancias

def resonancias(a, b, f, x0, ERROR):
    MAX_ITER = 10000
#    sol_newton = newton_raphson(f, x0, ERROR, MAX_ITER)
    sol_bisec = biseccion(f, a, b, ERROR, MAX_ITER)
    print("Solucion de la ecuacion no lineal:")
#    print("(newton)x = {}\tf({}) = {}".format(sol_newton, sol_newton, f(sol_newton)))
    print("(bisec)x = {}\tf({}) = {}".format(sol_bisec, sol_bisec, f(sol_bisec)))
    
    # Grafica de la funcion x0
    DISC = 1000
    x = np.linspace(a, b, DISC)
    plt.plot(x, f(x), label="f(x)")
#    plt.plot(sol_newton, f(sol_newton), "o", label="sol newton")
    plt.plot(sol_bisec, f(sol_bisec), "o", label="sol biseccion")
    plt.grid()
    plt.legend()
    plt.show()


    opt = input("Hallar otra raiz(yes/no): ")
    options = ["yes", "no"]
    opt.lower()
    while opt not in options:
        opt = input("(Opcion invalida) Hallar otra raiz(yes/no): ")
        opt.lower()
    
    if opt == "yes":
        frec_resonancias()
    else:
        sys.exit(1)

    

def inicializacion():

    print("Ingrese intervalo [a, b] a analizar")
    a0 = np.longdouble(input("Ingrese a(> 0): "))
    b0 = np.longdouble(input("Ingrese b: "))
    
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
    n = -(C2+C3)/(C1**2*C2**2*C3**2) + (6*C2*C3+4*C2**2+2*C3**2)*L3/(C1*C2**4*C3**2) - 2*(C2+C3)*R3**2/(C1*C2**3*C3**2) - R3**2/(C1*C2**4)
    l = lambda w: a(w)*(2*(C2+C3)*L2/(C1*C2**2*C3**2) + (2*C2+C3)*L3/C1**2*C2**2*C3 - R3**2/(C1**2*C2) - 2*R2*R3/(C1*C2**2)) \
                -(6*C2*C3+6*C2**2+C3**2)*L3**2/(C1*C2**4*C3**2) - R3**4/(C1*C2**2) + 2*(C2+C3)*L3*R3**2/(C1*C2**3*C3**2)
    p = lambda w: a(w)**2*(L2/C1**2 - R2**2/C1) + \
            ((2*L2*R3**2)/(C1*C2) -L3**2/(C1**2*C2) - 2*(2*C2+C3)*L2*L3/(C1*C2**2*C3))*a(w) +\
            2*(2*C2+C3)*L3**3/(C1*C2**3*C3) - 2*L3**2*R3**2/(C1*C2**2)

    r = lambda w: L1*a(w)**2*b(w) - L2**2*a(w)**2/C1 + 2*L2*L3**2*a(w)/(C1*C2) - L3**4/(C1*C2**2)

    f = lambda w: m + n*w**2 + l(w)*w**4 + p(w)*w**6 + r(w)*w**8

    """ Grafico para eleccion del punto de inicio del metodo de newton
    x = np.linspace(a0, b0, 1000)
    plt.plot(x, f(x), label="f(x)")
    plt.title("Identique la lozalización de la raices")
    plt.grid()
    plt.legend()
    plt.show()
"""
    
#    x0 = float(input("Punto de inicio(metodo de newton, SUG: valor ente {} y {}): ".format(a0, b0)))
    ERROR = float(input("Maxima cota de error: "))

    x0 = 0 # provisional (punto de inicio del metodo de newton
    return a0, b0, f, x0, ERROR
   
def frec_resonancias():
    a0, b0 , f , x0, ERROR = inicializacion()
    resonancias(a0, b0, f, x0, ERROR)

frec_resonancias()
