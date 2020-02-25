#Exercise 5.21: Electric Field of a charge distribution
#Written By: Nina Inman

import numpy as np
import matplotlib.pyplot as plt

import plotly.figure_factory as ff
import plotly.graph_objects as go

#go into interactive mode
plt.ion()

#Use Latex fonts
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

#a stepsize.... ???
h=0.01

#defines V from a distance r from the origin
def f(x, y, q):
    return q/(4 * np.pi * (8.85 * 10**(-12)) * (np.sqrt((x)**2 + (y)**2)))

#definitions of functions that calculate partial derivatives:
#sheet is two dimensions so we only need x and y...
def EField(x, y, q, h):
    Ex = (f((x + h/2), y, q) - f((x - h/2), y, q))/h
    Ey = ((f(x, (y + h/2), q) - f(x, (y - h/2), q))/h)
    #find the direction of the E field...
    mag = np.abs(Ex + Ey)
    return mag, Ex, Ey

####PART A####
#calculating the electric potential form two points on a plane
#initalizes the 2D array of potentials
potential = np.zeros([100, 100], float)
for x in range(100):
    for y in range(100):
        potential[x][y] = f(x - 55, y - 50, 1) + f(x - 45, y - 50, -1)

print(potential[0][0])
plt.hot()
plt.imshow(potential,  label = 'V')
plt.title('Electric Potential of Two Points')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
plt.savefig('exercise-521-potential.pdf',format='pdf')
plt.clf()


####PART B####
#initalizing a 2D arrays to save values from del...
field = np.zeros([100, 100], float)
x_dir = np.zeros([100, 100], float)
y_dir = np.zeros([100, 100], float)

run = 0
if (run == 1):
    #lets get ready to make a vector field:
    x_pos = []
    y_pos = []

    #calculating the total E field of 2 points 10cm apart...
    for x in range(100):
        for y in range(100):
            total_field1, Ex1, Ey1 = EField(x - 55, y - 50, 1, h)
            total_field2, Ex2, Ey2 = EField(x - 45, y - 50, -1, h)
            field[x][y] = total_field1 + total_field2
            x_dir[x][y] = Ex1 + Ex2
            y_dir[x][y] = Ey1 + Ey2


    plt.quiver(x_dir, y_dir, color = 'blue')
    plt.imshow(field,  origin = "lower")
    plt.title('Electric Field from Two Point Charges')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    plt.savefig('exercise-521-field.pdf',format='pdf')

####PART C####
def distribution(x, y):
    return 100 * cp.sin((2*np.pi*x)/0.1) * np.sin((2*np.pi*y))
#lets use guassian qudrature for the double integral...
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

def dist_V(a, b, N):
    #running it for 0 to 0.1(length of the distribution):
    x, w = gaussxwab(N, a, b)
    sum = 0
    step = 0.01

    #calculating the potential with the integral:
    for i in range(20):
        for j in range(20):
            sum = sum + (w[j] * w[i]) * disribution(x[i], x[j])

#Now to find the Electric field, varying r:
def EField(x, y, h):
    Ex = (distribution((x + h/2), y) - distribution((x - h/2), y))/h
    Ey = ((distribution(x, (y + h/2)) - distribution(x, (y - h/2)))/h)
    #find the direction of the E field...
    mag = np.abs(Ex + Ey)
    return mag, Ex, Ey
