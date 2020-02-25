#Exercise: 5.12 Stefan-Boltzmann constant
#Written By: Nina Inman

import numpy as np
import matplotlib.pyplot as plt

#go into interactive mode
plt.ion()

#Use Latex fonts
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

#initializing constants:
boltzmann = 1.38 * 10**(-23)
h_bar = 1.054 * 10**-34
c = 2.99 * 10**8

#Defining the function over which I'm integrating with a variable replacement:
# x = ((z)/(1 - z**2))
#dx = (1 + z**2)/(1 - z**2)**2
def f(z):
    if(z == 0):
        return 0
    elif(z == 1):
        return 0
    else:
        return ((((z)/(1 - z**2))**3)/(np.exp(((z)/(1 - z**2))) - 1))*((1 + z**2)/(1 - z**2)**2)

#defines a function for Simpson's Rules:
def Simpson(a, b, N):
    h = (b-a)/N
    sum = 0.
    x = 0.

    for i in range(1, N - 1, 1):
        #calculating the value of x for each iteration
        x = a + i * h
        #when the iteration is even, multiply by 2
        if(i % 2 == 0):
            sum = sum + (f(x))*2
        #when the iteration is odd, multiply by 4
        if(i % 2 == 1):
            sum = sum + (f(x))*4
    sum = (sum  + f(b))*(h/3)
    return sum, h

#integrating over infinity with a variable replacement:
a = 0.
b = 1.
N = 50

W, h = Simpson(a, b, N)
print("The value of the integral is given by:", W)

####PART C####
delta = W * ((boltzmann**4)/(4 * np.pi**2 * c**2 * h_bar**3))
print("The value of the Stefan-Boltzmann constant is", delta)
