import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
import multiprocessing as mp

def temporal_inc(l1,l2,t1,t2,dt):
    t = np.arange(t1,t2+dt,dt)
    return l1 + ((l2 - l1)*(t - t1)/(t2 - t1))

def HS(x,x0,n,l):
    return (1+l*((x/x0)**n))/(1+((x/x0)**n))

def integration_const(p,time,index,lamdas,l_temporal,n_periods):
    dt = time[1] - time[0]
    points = time.size - 1

    A = np.empty(points + 1)
    B = np.empty(points + 1)
    C = np.empty(points + 1)
    D = np.empty(points + 1)
    E = np.empty(points + 1)

    A[0] = p['Aic'][index]
    B[0] = p['Bic'][index]
    C[0] = p['Cic'][index]
    D[0] = p['Dic'][index]
    E[0] = p['Eic'][index]

    lDtoC = lamdas[0][0]
    lCtoD = lamdas[0][1]

    j = 0
    z = 0
    period = 0
    for i in range(1, points + 1):
        if (period < n_periods) and (100 + 300*period <= time[i] <= 150 + 300*period): #100 - 150 , 400 - 450 , 700 - 750
            lDtoC = l_temporal[0][j]
            lCtoD = l_temporal[1][j]
            j+=1
        if (period < n_periods) and (250 + 300*period <= time[i] <= 300 + 300*period): #250 - 300 , 550 - 600 ,
            lDtoC = l_temporal[0][j-1]
            lCtoD = l_temporal[1][j-1]
            j-=1
            if time[i] == 300*(1+period):
                period+=1
                j = 0

        A[i] = A[i-1] + dt * (p['gA'] * HS(C[i-1],p['C0'],p['nCtoA'],p['lCtoA']) - p['kA'] * A[i-1])
        B[i] = B[i-1] + dt * (p['gB'] * HS(A[i-1],p['A0'],p['nAtoB'],p['lAtoB']) - p['kB'] * B[i-1])
        C[i] = C[i-1] + dt * (p['gB'] * HS(B[i-1],p['B0'],p['nBtoC'],p['lBtoC']) * HS(D[i-1],p['D0'],p['nDtoC'],lDtoC) - p['kC'] * C[i-1])

        D[i] = D[i-1] + dt * (p['gD'] * HS(E[i-1],p['E0'],p['nEtoD'],p['lEtoD']) * HS(C[i-1],p['C0'],p['nCtoD'],lCtoD) - p['kD'] * D[i-1])
        E[i] = E[i-1] + dt * (p['gE'] * HS(D[i-1],p['D0'],p['nDtoE'],p['lDtoE']) - p['kE'] * E[i-1])

    return A,B,C,D,E


T = 1000
dt = 0.01
time = np.arange(0.0,T+dt,dt)

periods = 3

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
l_values = [[0.1,0.1],[0.9,0.9],[0.1,0.9],[0.9,0.1]]
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


if os.path.exists('C:\project\Summer_2022_DrMKJolly\RS_TS\size_20.npy'):
    arr = np.load('size_20.npy')
    p['Aic'] = arr[0]
    p['Bic'] = arr[1]
    p['Cic'] = arr[2]
    p['Dic'] = arr[3]
    p['Eic'] = arr[4]
else:
    #Initial Conditions
    p['Aic'] = np.random.random(size = 20) * (1.5 * (p['gA']/p['kA']))
    p['Bic'] = np.random.random(size = 20) * (1.5 * (p['gB']/p['kB']))
    p['Cic'] = np.random.random(size = 20) * (1.5 * (p['gC']/p['kC']))
    p['Dic'] = np.random.random(size = 20) * (1.5 * (p['gD']/p['kD']))
    p['Eic'] = np.random.random(size = 20) * (1.5 * (p['gE']/p['kE']))

    arr = np.array([p['Aic'],p['Bic'],p['Cic'],p['Dic'],p['Eic']])

    np.save('size_20.npy', arr)

combs = [[l_values[i],l_values[j]] for i in range(4) for j in range(4) if i!=j]
#t_change = np.arange(100.0,150+dt,dt)

for lamdas in combs:
    l_temporal = [temporal_inc(lamdas[0][0],lamdas[1][0],100,150,dt),temporal_inc(lamdas[0][1],lamdas[1][1],100,150,dt)]
    f, ax = plt.subplots(2,3,figsize = (16,12), sharex = True, sharey = True)
    for l in range(20):
        A,B,C,D,E = integration_const(p, time, l, lamdas, l_temporal, periods)
        ax[0,0].plot(time, A, linewidth = 1.5, label = str(p['Aic'][l]))
        ax[0,0].set_ylabel("A",fontsize=14)
        ax[0,1].plot(time, B, linewidth = 1.5, label = str(p['Bic'][l]))
        ax[0,1].set_ylabel("B",fontsize=14)
        ax[0,2].plot(time, C, linewidth = 1.5, label = str(p['Cic'][l]))
        ax[0,2].set_ylabel("C",fontsize=14)
        ax[1,0].plot(time, D, linewidth = 1.5, label = str(p['Dic'][l]))
        ax[1,0].set_ylabel("D",fontsize=14)
        ax[1,1].plot(time, E, linewidth = 1.5, label = str(p['Eic'][l]))
        ax[1,1].set_ylabel("E",fontsize=14)
    plt.tick_params(axis='y',labelsize=12,rotation=90)
    plt.tick_params(axis='x',labelsize=12)
    f.suptitle('lDtoC, lCtoD: ' + str(lamdas[0]) + '_' +str(lamdas[1]))
    ax[1,1].set_xlabel("Time (in weeks)",fontsize=14)
    f.tight_layout()
    plt.savefig('plot_'+'lDtoC_lCtoD_' + str(lamdas[0]) + '_' +str(lamdas[1])+'_50_'+str(periods)+'.jpeg',dpi=1000)
    plt.close()
