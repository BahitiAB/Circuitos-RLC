"""
Solucion de un sistema de odes(generado por un circuito rlc)
"""

import numpy as np


def runge_kutta_sist(x0, f, h, a, b):
    """
    Resuelve el sistema de ecuaciones diferenciales
    x' = f(t, x), x(t) = [x1(t), x2(t), ..., xm(t)], donde xr(t): R-> R^{p}, p natural
    input:
        x0: valor de inicio (manteniendo la estructura de x(t))
        f: funcion vectorial (camino)
        a, b: cota inferior y superior (resp) del intervalo [a, b]
        h: tamaño de paso
    output:
        x = [x0, ..., xk, ... xn]: xk= aproximacion de la solucion del sistema de odes tk = a + k*h.  
    """

    t = np.arange(a, b, h)
    x = np.array([x0])  # Solucion del sistema de odes

    for k, tk in enumerate(t):
        k1 = f(tk , x[k])
        k2 = f(tk + h/2, x[k]+h*k1/2)
        k3 = f(tk + h/2, x[k]+h*k2/2)
        k4 = f(tk + h, x[k]+h*k3)
        xx = x[k] + h*(k1+2*k2+2*k3+k4)/6 # xx = x_{k+1}
        x = np.append(x, xx) 
        
    return t, x


def f(t, x):
    z = x[0, :]
    q = x[1, :]

    # definiendo g(t, z, q) = [g1, g2, g3]
    
    V = lambda t: np.sin(t) # Voltaje variable
    L1 = 1.0    # Asignar valor real del inductancia
    L2 = 1.0    # Asignar valor real del inductancia
    L3 = 1.0    # Asignar valor real del inductancia

    C1 = 1.0    # Asignar valor real de la capacitancia
    C2 = 1.0    # Asignar valor real de la capacitancia
    C3 = 1.0    # Asignar valor real de la capacitancia


    R1 = 1.0    # Asignar valor real de la resistencia
    R2 = 1.0    # Asignar valor real de la resistencia
    R3 = 1.0    # Asignar valor real de la resistencia

    #NOTA: qk = q[k-1]
    
    
    g1 = lambda t, z, q: V(t)/L1 - ((L1+L2)/(C1*L1*L2))*q[0] + 1/(L2*C2)*q[1] - R1/L1*z[0] - (L2/R2-L1/R1)*(z[1]+z[2])
    g2 = lambda t, z, q: q[0]/(L2*C1) - (L2+L3)/(C2*L2*L3)*q[1] + 1/(C3*L3)*q[2] - (R2/L2)*(z[1]+z[2])
    g3 = lambda t, z, q: q[1]/(L3*C2) - q[2]/(C3*L3) - (R3/L3)*z[2] 
    
    g = lambda t, z, q: np.array([g1(t, z, q), g2(t, z, q), g3(t, z, q)])

    f1 = lambda t, z, q: g(t, z, q)
    f2 = lambda t, z, q: z

    f = lambda t, z, q: np.array([f1(t, z, q), f2(t, z, q)])

    return f(t, z, q)




#NOTA: x0 = [z0, q0] = [q'(0), q(0)] = [i(0), q(0)], donde i = (i1, i2, i3) y q = (q1, q2, q3)
x0 = np.array([[0, 0, 0], [0, 0, 0]])

# [a, b] es el intervalo de tiempo en el cual se estara analizando el circuito
a = 0   # Inicializar valor real
b = 1   # Inicializar valor real
h = 0.01    # tamaño de paso en el metodo de runge-kutta

x = runge_kutta_sist(x0, f, h, a, b)