#Implementation of Jon's Example_plot.py but using the new OOP code
from Fncts import *

import numpy as np
from scipy import linalg
import math

#packages for plotting
import matplotlib
import matplotlib.pyplot as plt


V0 = 0.1
reservoir_rate = 10**7 #This is the default rate constant for electron flow INTO the reservoirs. This can be changed directly using the last argument of the Add_reservoir() function


slope = 0.15 #this is the "slope" of the energy landscapes (i.e. the difference in reduction potentials of neighboring cofactors)

cofactor_distance = 10 #This is the distance between neighboring cofactors


N = 100 #The number of points to be plotted
res2emin = -0.1 #The range of energies of the 2-electron (D) reservoir to be plotted over
res2emax = 0.1
dx = (res2emax-res2emin)/N #energy step size


data1  = [] #initiate arrays to store data
data2 = []
data3 = []


for n in range(N):
    Net=Network()

    #Initialize the 4 cofactors
    #1st element in list is the reduction potential from 0e- to 1e-
    D = Cofactor("D", [-0.4, 0.4])
    B = Cofactor("B", [-0.4, 0.4])
    H1 = Cofactor("H1", [0.4 - slope*1])
    H2 = Cofactor("H2", [0.4  - slope*2])
    L1 = Cofactor("L1", [-0.4 + slope*1])
    L2 = Cofactor("L2", [-0.4 + slope*2])



    Net.addCofactor(D)
    Net.addCofactor(B)
    Net.addCofactor(H1)
    Net.addCofactor(H2)
    Net.addCofactor(L1)
    Net.addCofactor(L2)

    #Connect the fully reduced form of the bifurcating cofactor to the proximal acceptors

    Net.addConnection(D, B, 7)
    Net.addConnection(H1, H2, 10)
    Net.addConnection(B, H1, 10)
    Net.addConnection(B, L1, 10) #this is a short circuit
    Net.addConnection(B, H2, 20)
    Net.addConnection(B, L2, 20) #this is a short circuit

    #Connect cofactors further down the branches

    Net.addConnection(L1, L2, 10)

    #Add short-circuit rates directly between branches
    Net.addConnection(H1, L1, 20) #this is a short circuit
    Net.addConnection(H2, L1, 30) #this is a short circuit
    Net.addConnection(H1, L2, 30) #this is a short circuit
    Net.addConnection(H2, L2, 40) #this is a short circuit
    
    #Add rate constants between D and high- and low-potential branches
    Net.addConnection(D, H1, math.sqrt(10**2 + 7 **2))
    Net.addConnection(D, H2, math.sqrt(20**2 + 7 **2))
    Net.addConnection(D, L1, math.sqrt(10**2 + 7 **2))
    Net.addConnection(D, L2, math.sqrt(20**2 + 7 **2))

    #Connect bifurcating and terminal cofactors to reservoirs
    Net.addReservoir("Two-Electron Reservoir", D, 2, 2,res2emin + dx*n, reservoir_rate)
    Net.addReservoir("High Potential Reservoir", H2, 1, 1, 0, reservoir_rate)
    Net.addReservoir("Low Potential Reservoir",L2, 1, 1, 0, reservoir_rate)

    #Construct Adjacency matrix
    Net.constructAdjacencyMatrix()

    #Construct Rate matrix
    Net.constructRateMatrix()

    #Add the concerted transfer from D to B
    kc = 10**3 #concerted two-electron rate constant
    #Net.addMultiElectronConnection(D, B, 2, 0, 2, kc)

    #Set the initial population
    pop_init=np.zeros(Net.num_state)
    pop_init[0]=1
    t = 20

    pop = Net.evolve(t, pop_init)


    #def getReservoirFlux(self, reservoir_id: int, cofactor: Cofactor, initial_redox: int, final_redox: int, pop: np.array) -> float:
    #Plot the reservoir flux over time
    data1.append(Net.getReservoirFlux("Two-Electron Reservoir", pop))
    data2.append(Net.getReservoirFlux("High Potential Reservoir", pop))
    data3.append(Net.getReservoirFlux("Low Potential Reservoir", pop))


x = np.linspace(res2emin*1000, res2emax*1000, N) #*1000 to convert to meV

plt.rc('font', family='DejaVu Sans')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
plt.rc('text')

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)
ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
plt.plot(x, data1, '#91009B', linestyle='--')
plt.plot(x, data2,'#2D00C8')
plt.plot(x, data3,'#B40005')
ax.set_xlabel('$\Delta G_{bifurc}$ (meV)',size='x-large')
ax.set_ylabel('Flux (Sec$^{-1}$)',size='x-large')
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(left=0.2)

plt.savefig("3D.svg")
