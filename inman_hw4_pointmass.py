#Exercise: 5.14 Gravitational Force from a uniform sheet
#Written By: Nina Inman

import numpy as np
import matplotlib.pyplot as plt

#go into interactive mode
plt.ion()

#Use Latex fonts
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

#constants
G = 6.67 * 10**(-11)

#defining the function over which we are integrating:
def density(x, y, z):
    return 1/((x**2 + y**2 + z**2)**(3/2))
#vary over x and y with list
#definition of Guassian Quadrature (Source: Newmann):
def gaussxw(N):
    #approximation to roots of Legendre polynomial:
    a = np.linspace(3, 4*N-1, N)/(4*N+2)
    x = np.cos(np.pi * a+1 / (8*N*N*np.tan(a)))

    #finding roots using Newton's method
    epsilon = 1 * (10**-15)
    delta = 1.
    while delta > epsilon:
        p0 = np.ones(N, float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(np.abs(dx))

    #calculate the weights
    w = 2 * (N + 1) * (N + 1)/(N * N *(1 - x*x)*dp*dp)

    return x, w

def gaussxwab(N, a, b):
    x, w = gaussxw(N)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

#initializing variables
#limits of integration for y:
a = -5
b = 5
N = 50
z = 0
forcy = []
z_vals = []

#first we can integrate from 0 to 1 for F(y)
#integrates through each possible point in the block with the weights found in guassian quadrature
x, w = gaussxwab(N, a, b)
step = 10/100
for k in range(0, 100, 1):
    force = 0
    #adds up the force from each point on the sheet for a certain z
    for i in range(0, N, 1):
        for j in range(0, N, 1):
            force = force + G * z * (10000/100) *(w[i] * w[j] * density(x[i], x[j], z))
    forcy.append(force)
    z_vals.append(z)
    z += step

plt.plot(z_vals, forcy, label = "$F_{z}$")
plt.xlabel('Z (distance of point mass from sheet) (m)')
plt.ylabel('Gravitational Force (N)')
plt.title('Force as a Function of Time')
plt.legend(loc = 1)
plt.show()
plt.savefig('forceofsheet.pdf',format='pdf')
