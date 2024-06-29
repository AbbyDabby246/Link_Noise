import numpy as np
import os
import random
from itertools import permutations
import matplotlib.pyplot as plt
#import multiprocessing as mp

no_of_cpu_for_use=5

def HS(x,x0,n,l):
    return (1+l*((x/x0)**n))/(1+((x/x0)**n))

def step(x, tau):
    return 1 * (x > tau)

def integration_const(p,time,time2,index,p_index,noise, lamda):
    dt = time[1] - time[0]
    points = time.size - 1
    dt2 = time2[1] - time2[0]

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

    lEtoD = p['lEtoD'][p_index]
    lDtoE = p['lDtoE'][p_index]
    lAtoB = p['lAtoB'][p_index]
    lBtoC = p['lBtoC'][p_index]
    lCtoA = p['lCtoA'][p_index]
    
    lDtoC = lamda[0]
    lCtoD = lamda[1]

    for i in range(1, points + 1):
        if time[i] % 1000 == 0:
            print(time[i])

        
        if time[i] in time2:
            lEtoD = max(0,min(1,(lEtoD + random.gauss(0,noise))))#max(lDtoC + random.uniform(0.1,0.5),0) if random.random()<0.5 else max(lDtoC - random.uniform(0.1,0.5),0)
            lDtoE = max(0,min(1,(lDtoE + random.gauss(0,noise))))
            lAtoB = max(0,min(1,(lAtoB + random.gauss(0,noise))))
            lBtoC = max(0,min(1,(lBtoC + random.gauss(0,noise))))
            lCtoA = max(0,min(1,(lCtoA + random.gauss(0,noise))))#max(lCtoD + random.uniform(0.1,0.5),0) if random.random()<0.5 else max(lCtoD - random.uniform(0.1,0.5),0)
            #print(time[i], lEtoD, lDtoE,lAtoB,lBtoC,lCtoA)
        #lamda_random.append([lDtoE,lDtoE,lAtoB,lBtoC,lCtoA])


        A[i] = A[i-1] + dt * (p['gA'][p_index] * HS(C[int(i - 1 - p['tau'] * step(i, p['tau']))], p['C0'][p_index], p['nCtoA'][p_index], lCtoA) - p['kA'][p_index] * A[i-1])
        B[i] = B[i-1] + dt * (p['gB'][p_index] * HS(A[int(i - 1 - p['tau'] * step(i, p['tau']))], p['A0'][p_index], p['nAtoB'][p_index], lAtoB) - p['kB'][p_index] * B[i-1])
        C[i] = C[i-1] + dt * (p['gB'][p_index] * HS(B[int(i - 1 - p['tau'] * step(i, p['tau']))], p['B0'][p_index], p['nBtoC'][p_index], lBtoC) * HS(D[i-1],p['D0'][p_index],p['nDtoC'],lDtoC) - p['kC'][p_index] * C[i-1])

        D[i] = D[i-1] + dt * (p['gD'][p_index] * HS(E[int(i - 1 - 0 * step(i, 0))], p['E0'][p_index], p['nEtoD'][p_index], lEtoD) * HS(C[i-1],p['C0'][p_index],p['nCtoD'],lCtoD) - p['kD'][p_index] * D[i-1])
        E[i] = E[i-1] + dt * (p['gE'][p_index] * HS(D[int(i - 1 - 0 * step(i, 0))], p['D0'][p_index], p['nDtoE'][p_index], lDtoE) - p['kE'][p_index] * E[i-1])

    return A,B,C,D,E



T = 10000
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
dt2 = 0.05
time2 = np.arange(1000.0,T+dt,dt2).round(2)
noises = [0.001]#[0.0001,0.001,0.002,0.005] #variance of normal dist
n_points = time.size
print(time2)

folder = "C:/project/Summer_2022_DrMKJolly/RS_TS/plots/dynamic/random_indi/delay/"

p = {}

#production paramters
p['gA'] = [20.493077, 26.721114, 55.061022, 39.731889, 74.811574]
p['gB'] = [22.988018, 84.273707, 33.813016, 85.411491, 92.607457]
p['gC'] = [40.181269, 80.608903, 55.735976, 29.439667, 26.317553]
p['gD'] = [91.947511, 50.212433, 64.102184, 83.975461, 39.691404]
p['gE'] = [38.696894, 95.256533, 25.889225, 66.912087, 77.538958]

#Threshold constants
p['A0'] = [ 8.766943, 44.061648, 22.812194,  3.313576, 51.741789]
p['B0'] = [28.432453, 29.876363, 44.228931, 21.493818, 23.456629]
p['C0'] = [16.655017, 42.603391, 20.022012,  4.426666, 30.799155]
p['D0'] = [53.920926, 14.888436, 31.018099, 53.370193,  3.534087]
p['E0'] = [36.738347, 35.491987,  7.507220, 10.120514, 42.641427]

#Degradation constants
p['kA'] = [0.127069, 0.111813, 0.387255, 0.546236, 0.433251]
p['kB'] = [0.349802, 0.781383, 0.172455, 0.483953, 0.615320]
p['kC'] = [0.287023, 0.638218, 0.522554, 0.749755, 0.122211]
p['kD'] = [0.699470, 0.390561, 0.829410, 0.207175, 0.881824]
p['kE'] = [0.442699, 0.444709, 0.801409, 0.498639, 0.431325]

#Fold Change
l_values = [0.1, 0.9, 0.5]#[0.1,0.1],[0.9,0.9],[0.1,0.9],[0.9,0.1]]
p['lAtoB'] = [0.025056, 0.057252, 0.012882, 0.019468, 0.021675]
p['lBtoC'] = [0.012458, 0.010081, 0.012339, 0.015197, 0.011492]
p['lCtoA'] = [0.010947, 0.012998, 0.050211, 0.029954, 0.013929]
p['lDtoE'] = [0.011075, 0.028405, 0.017531, 0.015665, 0.011402]
p['lEtoD'] = [0.021791, 0.018186, 0.018946, 0.020940, 0.018222]


#Hill Coeff
p['nAtoB'] = [2.000000, 6.000000, 3.000000, 2.000000, 3.000000]
p['nBtoC'] = [4.000000, 4.000000, 5.000000, 5.000000, 3.000000]
p['nCtoA'] = [6.000000, 4.000000, 4.000000, 5.000000, 5.000000]
p['nDtoE'] = [3.000000, 6.000000, 6.000000, 2.000000, 5.000000]
p['nEtoD'] = [5.000000, 3.000000, 3.000000, 5.000000, 6.000000]
p['nDtoC'] = 5.000000
p['nCtoD'] = 5.000000

p['tau'] = 0.0

n = 50

lamdas = list(permutations(l_values, r = 2))
lamdas.append((0.1, 0.1))
lamdas.append((0.9, 0.9))
lamdas.append((0.5, 0.5))

for i in range(0, 1, 1):

    if os.path.exists(folder + 'npy_files/size_50_random'+str(i)+'.npy'):
        arr = np.load(folder + 'npy_files/size_50_random'+str(i)+'.npy')
        p['Aic'] = arr[0]
        p['Bic'] = arr[1]
        p['Cic'] = arr[2]
        p['Dic'] = arr[3]
        p['Eic'] = arr[4]
    else:
        #Initial Conditions
        p['Aic'] = np.random.random(size = n) * (1.5 * (p['gA'][i]/p['kA'][i]))
        p['Bic'] = np.random.random(size = n) * (1.5 * (p['gB'][i]/p['kB'][i]))
        p['Cic'] = np.random.random(size = n) * (1.5 * (p['gC'][i]/p['kC'][i]))
        p['Dic'] = np.random.random(size = n) * (1.5 * (p['gD'][i]/p['kD'][i]))
        p['Eic'] = np.random.random(size = n) * (1.5 * (p['gE'][i]/p['kE'][i]))

        arr = np.array([p['Aic'],p['Bic'],p['Cic'],p['Dic'],p['Eic']])

        np.save(folder + 'npy_files/size_50_random'+str(i)+'.npy', arr)

    for noise in noises:
        #arr_dynamics = []
        print(noise)
        
        #lamdas_random = [random.choice(lamdas)]
        
        for lamda in [(0.1,0.9)]:
            
            arr_lamda = []
            
            f, ax = plt.subplots(2,3,figsize = (16,12), sharex = True, sharey = True)
            
            for l in range(n):
                A,B,C,D,E = integration_const(p, time, time2, l, i, noise, lamda)
                arr_lamda.append(np.array([A,B,C,D,E]))
                
                ax[0,0].plot(time, A, linewidth = 1.5, label = str(p['Aic'][l]), alpha = 0.3)
                ax[0,0].set_ylabel("A",fontsize=14)
                ax[0,1].plot(time, B, linewidth = 1.5, label = str(p['Bic'][l]), alpha = 0.3)
                ax[0,1].set_ylabel("B",fontsize=14)
                ax[0,2].plot(time, C, linewidth = 1.5, label = str(p['Cic'][l]), alpha = 0.3)
                ax[0,2].set_ylabel("C",fontsize=14)
                ax[1,0].plot(time, D, linewidth = 1.5, label = str(p['Dic'][l]), alpha = 0.3)
                ax[1,0].set_ylabel("D",fontsize=14)
                ax[1,1].plot(time, E, linewidth = 1.5, label = str(p['Eic'][l]), alpha = 0.3)
                ax[1,1].set_ylabel("E",fontsize=14)
            
            plt.tick_params(axis='y',labelsize=12,rotation=90)
            plt.tick_params(axis='x',labelsize=12)
            f.suptitle('lDtoC, lCtoD: ' + str(lamda[0]) + '_' +str(lamda[1]))
            ax[1,1].set_xlabel("Time (in weeks)",fontsize=14)
            f.tight_layout()
            plt.savefig(folder + 'plot_delayDynamicsCoupled_noise_' + str(noise) + '_dt2_' + str(dt2) + '_time_' + str(T) + '_param_' + str(i) + '_delay_' + str(p['tau']) + '_lDtoC_lCtoD_' + str(lamda[0]) + '_' + str(lamda[1]) +'.jpeg',dpi=1000)
            plt.show()
            plt.close()
            
            
        
            #arr_dynamics.append(arr_lamda)
            
        file = 'value_storage/delayDynamicsCoupled_noise_' + str(noise) + '_dt2_' + str(dt2) + '_time_' + str(T) + '_param_' + str(i) + '_delay_' + str(p['tau']) + '.npy'
        #np.save(folder + file, np.array(arr_dynamics))
