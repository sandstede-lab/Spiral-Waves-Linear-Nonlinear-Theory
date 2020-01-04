import math
import numpy as np

file_input = 'Barkley_transport_ic_original.txt';
file_output = 'Barkley_transport_ic_perturbed.txt';

with open(file_input) as myfile:
    header = [next(myfile) for x in range(5)]

U = np.loadtxt(open(file_input, "rb"), usecols=0, skiprows=10)
V = np.loadtxt(open(file_input, "rb"), usecols=1, skiprows=10)
N = int(math.sqrt(U.size))
V = np.reshape(V, (N, N))

L = 100
epsilon = 0.1
distance = 20
spread = 3

y = np.linspace(-L,L,N)
[X,Y] = np.meshgrid(y,y)
R = np.sqrt(X**2 + Y**2)
P = epsilon/np.cosh((R - distance)/spread)
VP = V+P

U = np.reshape(U, (N*N, 1))
V = np.reshape(VP, (N*N, 1))
Z = np.column_stack((U, V))

with open(file_output, 'w+') as filehandle:
    for listitem in header:
        filehandle.write('%s' % listitem)
    np.savetxt(filehandle, Z)
