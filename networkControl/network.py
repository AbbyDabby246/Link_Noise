import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
import multiprocessing as mp

def HS(x,x0,n,l):
    return (1+l*((x/x0)**n))/(1+((x/x0)**n))

def integration_const(p,time,z):
    dt = time[1] - time[0]
    points = time.size - 1

    A = np.empty(points + 1)
    B = np.empty(points + 1)
    C = np.empty(points + 1)
    D = np.empty(points + 1)
    E = np.empty(points + 1)

    A[0] = p['Aic'][z]
    B[0] = p['Bic'][z]
    C[0] = p['Cic'][z]
    D[0] = p['Dic'][z]
    E[0] = p['Eic'][z]

    for i in range(1, points + 1):
        A[i] = A[i-1] + dt * (p['gA'] * HS(C[i-1],p['C0'],p['nCtoA'],p['lCtoA']) - p['kA'] * A[i-1])
        B[i] = B[i-1] + dt * (p['gB'] * HS(A[i-1],p['A0'],p['nAtoB'],p['lAtoB']) - p['kB'] * B[i-1])
        C[i] = C[i-1] + dt * (p['gB'] * HS(B[i-1],p['B0'],p['nBtoC'],p['lBtoC']) * HS(D[i-1],p['D0'],p['nDtoC'],p['lDtoC']) - p['kC'] * C[i-1])

        D[i] = D[i-1] + dt * (p['gD'] * HS(E[i-1],p['E0'],p['nEtoD'],p['lEtoD']) * HS(C[i-1],p['C0'],p['nCtoD'],p['lCtoD']) - p['kD'] * D[i-1])
        E[i] = E[i-1] + dt * (p['gE'] * HS(D[i-1],p['D0'],p['nDtoE'],p['lDtoE']) - p['kE'] * E[i-1])

    return A,B,C,D,E


T = 100
dt = 0.01
time = np.arange(0.0,T+dt,dt)

p = {}

#production paramters
p['gA'] = 20.493077
p['gB'] = 22.988018
p['gC'] = 40.181269
p['gD'] = 91.947511
p['gE'] = 38.696894

#Threshold constants
p['A0'] = 8.766943
p['B0'] = 28.432453
p['C0'] = 16.655017
p['D0'] = 53.920926
p['E0'] = 36.738347

#Degradation constants
p['kA'] = 0.127069
p['kB'] = 0.349802
p['kC'] = 0.287023
p['kD'] = 0.699470
p['kE'] = 0.442699

#Fold Change
l_values = [0.1,0.9,2]
p['lAtoB'] = 0.025056
p['lBtoC'] = 0.012458
p['lCtoA'] = 0.010947
p['lDtoE'] = 0.011075
p['lEtoD'] = 0.021791
#p['lDtoC'] = [0.1, 0.5, 0.9]
#p['lCtoD'] = [0.1, 0.5, 0.9]

#Hill Coeff
p['nAtoB'] = 2.000000
p['nBtoC'] = 4.000000
p['nCtoA'] = 6.000000
p['nDtoE'] = 3.000000
p['nEtoD'] = 5.000000
p['nDtoC'] = 5.000000
p['nCtoD'] = 5.000000

#Initial Conditions
p['Aic'] = np.random.uniform(size = 50) * (1.5 * (p['gA']/p['kA']))
p['Bic'] = np.random.uniform(size = 50) * (1.5 * (p['gB']/p['kB']))
p['Cic'] = np.random.uniform(size = 50) * (1.5 * (p['gC']/p['kC']))
p['Dic'] = np.random.uniform(size = 50) * (1.5 * (p['gD']/p['kD']))
p['Eic'] = np.random.uniform(size = 50) * (1.5 * (p['gE']/p['kE']))


for i in range(1):
    p['lDtoC'] = 0.8#l_values[i]
    for j in range(1):
        p['lCtoD'] = 0.7#l_values[j]
        f, ax = plt.subplots(2,3,figsize = (16,12), sharex = True, sharey = True)
        for z in range(50):
            A,B,C,D,E = integration_const(p, time, z)
            ax[0,0].plot(time, A, linewidth = 1.5, label = str(p['Aic'][z]))
            ax[0,0].set_ylabel("A",fontsize=14)
            ax[0,1].plot(time, B, linewidth = 1.5, label = str(p['Bic'][z]))
            ax[0,1].set_ylabel("B",fontsize=14)
            ax[0,2].plot(time, C, linewidth = 1.5, label = str(p['Cic'][z]))
            ax[0,2].set_ylabel("C",fontsize=14)
            ax[1,0].plot(time, D, linewidth = 1.5, label = str(p['Dic'][z]))
            ax[1,0].set_ylabel("D",fontsize=14)
            ax[1,1].plot(time, E, linewidth = 1.5, label = str(p['Eic'][z]))
            ax[1,1].set_ylabel("E",fontsize=14)
        plt.tick_params(axis='y',labelsize=12,rotation=90)
        plt.tick_params(axis='x',labelsize=12)
        f.suptitle('lDtoC: ' + str(p['lDtoC']) + ', lCtoD: ' + str(p['lCtoD']))
        ax[1,1].set_xlabel("Time (in weeks)",fontsize=14)
        f.tight_layout()
        plt.savefig('plot_'+str(p['lDtoC'])+'_'+str(p['lCtoD'])+'.jpeg',dpi=1000)
        plt.close()

'''
p['lDtoC'] = 1
p['lCtoD'] = 1
f, ax = plt.subplots(2,3,figsize = (16,12), sharex = True, sharey = True)
for z in range(5):
    A,B,C,D,E = integration_const(p, time, z)
    ax[0,0].plot(time/(24*7), A, linewidth = 1.5, label = str(p['Aic'][z]))
    ax[0,0].set_ylabel("A",fontsize=14)
    ax[0,1].plot(time/(24*7), B, linewidth = 1.5, label = str(p['Bic'][z]))
    ax[0,1].set_ylabel("B",fontsize=14)
    ax[0,2].plot(time/(24*7), C, linewidth = 1.5, label = str(p['Cic'][z]))
    ax[0,2].set_ylabel("C",fontsize=14)
    ax[1,0].plot(time/(24*7), D, linewidth = 1.5, label = str(p['Dic'][z]))
    ax[1,0].set_ylabel("D",fontsize=14)
    ax[1,1].plot(time/(24*7), E, linewidth = 1.5, label = str(p['Eic'][z]))
    ax[1,1].set_ylabel("E",fontsize=14)
plt.tick_params(axis='y',labelsize=12,rotation=90)
plt.tick_params(axis='x',labelsize=12)
f.suptitle('lDtoC: ' + str(p['lDtoC']) + ', lCtoD: ' + str(p['lCtoD']))
ax[1,1].set_xlabel("Time (in weeks)",fontsize=14)
f.tight_layout()
plt.savefig('plot_'+str(p['lDtoC'])+'_'+str(p['lCtoD'])+'.jpeg',dpi=1000)
plt.close()

'''


#A3,B3,C3,D3,E3 = integration_const(p, time,0,200,0,0,200,0)
#A2,B2,C2,D2,E2 = integration_const(p, time,0,0,200,200,0,0)
#A3,B3,C3,D3,E3 = integration_const(p, time,0,0,200,100,0,0)
#A4,B4,C4,D4,E4 = integration_const(p, time,0,0,0,0,0,0)
#A5,B5,C5,D5,E5 = integration_const(p, time,0,0,100,100,0,0)
'''
f, ax = plt.subplots(figsize=(4,3))
ax.plot(time,C1,'g-', linewidth = 2, label = 'C1')
ax.plot(time,C2,'b-', linewidth = 2, label = 'C2')
#ax.plot(time,C3,'r-', linewidth = 2, label = 'C3')
#ax.plot(time,C4,'k-', linewidth = 2, label = 'C4')
#ax.plot(time,C5,'y-', linewidth = 2, label = 'C5')

ax.tick_params(axis='y',labelsize=12,rotation=90)
ax.tick_params(axis='x',labelsize=12)
ax.set_ylabel("C and D",fontsize=14)
ax.set_xlabel("Time (in weeks)",fontsize=14)

plt.legend(loc = 'upper right', ncol=2)
plt.subplots_adjust(top=0.95, bottom=0.20, left=0.20, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig('C_'+str(p['lDtoC'][0])+str(p['lCtoD'][0])+'_2.jpeg',dpi=2000)
plt.close()

f, ax = plt.subplots(figsize=(4,3))
ax.plot(time,D1,'g-', linewidth = 2, label = 'D1')
ax.plot(time,D2,'b-', linewidth = 2, label = 'D2')
#ax.plot(time,D3,'r-', linewidth = 2, label = 'D3')
#ax.plot(time,D4,'k-', linewidth = 2, label = 'D4')
#ax.plot(time,D5,'y-', linewidth = 2, label = 'D5')

ax.tick_params(axis='y',labelsize=12,rotation=90)
ax.tick_params(axis='x',labelsize=12)
ax.set_ylabel("C and D",fontsize=14)
ax.set_xlabel("Time (in weeks)",fontsize=14)

plt.legend(loc = 'upper right', ncol=2)
plt.subplots_adjust(top=0.95, bottom=0.20, left=0.20, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig('D_'+str(p['lDtoC'][0])+str(p['lCtoD'][0])+'_2.jpeg',dpi=2000)
plt.close()
'''
