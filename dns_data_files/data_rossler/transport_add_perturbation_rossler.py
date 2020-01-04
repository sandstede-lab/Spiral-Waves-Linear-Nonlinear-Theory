import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import numpy as np

file_input = 'Rossler_transport_ic_original.txt';
file_output = 'Rossler_transport_ic_perturbed.txt';

with open(file_input) as myfile:
    header = [next(myfile) for x in range(5)]
    
U = np.loadtxt(open(file_input, "rb"), usecols=0, skiprows=5)
V = np.loadtxt(open(file_input, "rb"), usecols=1, skiprows=5)
W = np.loadtxt(open(file_input, "rb"), usecols=2, skiprows=5)
N = int(math.sqrt(U.size))
U = np.reshape(U, (N, N))

L = 300
epsilon = 0.2
distance = 100
spread = 50

y = np.linspace(-L,L,N)
[X,Y] = np.meshgrid(y,y)
R = np.sqrt(X**2 + Y**2)
P = epsilon/np.cosh((R - distance)/spread)
UP = U+P

U = np.reshape(UP, (N*N, 1))
V = np.reshape(V, (N*N, 1))
W = np.reshape(W, (N*N, 1))
Z = np.column_stack((U, V, W))

with open(file_output, 'w+') as filehandle:
    for listitem in header:
        filehandle.write('%s' % listitem)
    np.savetxt(filehandle, Z)
