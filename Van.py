#! /usr/bin/env python

'''Script que realiza RK3 para integrar la ecuación de van der Pol
y grafica el resultado'''

import numpy as np
import matplotlib.pyplot as plt

mu=1.844 #'''parametro de la ecuacion'''

def f (y,m):
    '''funcion a la que se le aplica RK3'''
    return (m,-y-mu*(y**2-1)*m)

def get_k1(y_n,m_n,h,f):
    '''calculo de k1'''
    f_n = f(y_n,m_n)
    return h*f_n[0],h*f_n[1]

def get_k2(y_n,m_n,h,f):
    '''calculo de k2'''
    k1 = get_k1(y_n,m_n,h,f)
    f_n = f(y_n+k1[0]/2.,m_n+k1[1]/2.)
    return k1,(h*f_n[0],h*f_n[1])

def get_k3(y_n,m_n,h,f):
    '''calculo de k3'''
    k1,k2 = get_k2(y_n,m_n,h,f)
    f_n = f(y_n-k1[0]+2*k2[0],m_n-k1[1]+2*k2[1])
    return k1,k2,(h*f_n[0],h*f_n[1])

def avanzar_rk3 (y_n,m_n,h,f):
    '''recibe los valores en el paso n-esimo de "y" y "m"
    y retorna los valores en el paso siguiente'''
    k1,k2,k3 = get_k3(y_n,m_n,h,f)
    y_n1 = y_n + 1./6. * (k1[0]+4*k2[0]+k3[0])
    m_n1 = m_n + 1./6. * (k1[1]+4*k2[1]+k3[1])
    return y_n1,m_n1
'''-------------------------------------------------------------------------------------------'''

'''condiciones iniciales m0=0 y0=0.1'''
m0 = 0
y0 = 0.1

n_pasos = 1000
h = 20*np.pi / n_pasos
y_1 = np.zeros(n_pasos)
m_1 = np.zeros(n_pasos)

y_1[0] = y0
m_1[0] = m0

for i in range(1,n_pasos):
    (y_1[i],m_1[i]) = avanzar_rk3(y_1[i-1],m_1[i-1],h,f)

plt.figure(1)
plt.clf
plt.plot(y_1,m_1,color="r",label="condiciones iniciales: dy/ds=0, y=0.1")
plt.xlabel('$y$', fontsize=20)
plt.ylabel("$\\frac{dy}{ds}$",fontsize=20)
plt.legend(loc='lower right',prop={'size':10})
plt.title("Oscilador de Van der Pol")
plt.savefig("Van der pol1.png")
'''condiciones iniciales m0=0 y0=4'''
m0 = 0
y0 = 4

n_pasos = 1000
h = 20*np.pi / n_pasos
y_2 = np.zeros(n_pasos)
m_2 = np.zeros(n_pasos)

y_2[0] = y0
m_2[0] = m0

for i in range(1,n_pasos):
    (y_2[i],m_2[i]) = avanzar_rk3(y_2[i-1],m_2[i-1],h,f)

plt.figure(2)
plt.clf
plt.plot(y_2,m_2,color="g",label="condiciones iniciales: dy/ds=0, y=4")
plt.xlabel('$y$', fontsize=20)
plt.ylabel("$\\frac{dy}{ds}$",fontsize=20)
plt.legend(loc='lower right',prop={'size':10})
plt.title("Oscilador de Van der Pol")
plt.savefig("Van der pol2.png")

plt.figure(3)
plt.clf
plt.plot(y_2,m_2,color="g",label="condiciones iniciales: dy/ds=0, y=4")
plt.plot(y_1,m_1,color="r",label="condiciones iniciales: dy/ds=0, y=4")
plt.xlabel('$y$', fontsize=20)
plt.ylabel("$\\frac{dy}{ds}$",fontsize=20)
plt.legend(loc='lower right',prop={'size':10})
plt.title("Oscilador de Van der Pol")
plt.savefig("Van der pol3.png")

plt.show()
