#! /usr/bin/env python

'''Script que realiza RK4 para integrar el sistema de Lorenz'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode

sigma = 10
beta = 8/3
rho = 28

def f(t,(x,y,z)):
    '''funcion que se integrara'''
    return [sigma*(y-x),x*(rho-z)-y,x*y-beta*z]

'''condiciones iniciales'''
x0 = 1
y0 = 1
z0 = 1

solucion=ode(f).set_integrator("dopri5")
solucion.set_initial_value([x0,y0,z0],0)

t1 = 100
dt = 0.01
pasos = t1/dt +1
i = 0
t_n = np.zeros(pasos)
sol_n = [np.zeros(pasos),np.zeros(pasos),np.zeros(pasos)]
while solucion.successful() and solucion.t < t1:
    t_n[i],(sol_n[0][i],sol_n[1][i],sol_n[2][i])=[solucion.t,(solucion.integrate(solucion.t+dt))]
    i += 1

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

ax.plot(sol_n[0],sol_n[1] ,sol_n[2] )

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title("Atractor de Lorenz")
fig.savefig("Lorenz.png")
fig.show()
