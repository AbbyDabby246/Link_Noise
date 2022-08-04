import numpy as np
import os
import random
import matplotlib.pyplot as plt
from statistics import stdev, mean
#import multiprocessing as mp

no_of_cpu_for_use=5

def check_thresh(l,thr,temp = 0):
    if temp:
        c = [n > thr for n in l[50000:]]
    else:
        c = [n < thr for n in l[50000:]]
    
    if sum(c) > 0:
        return c.index(1)
    return -1

def find_states(arr, a, n):
    unique = []
    for i in range(n):
        unique.append(arr[i][a][40000])
    unique = list(set(unique))
    unique_sorted = sorted(unique)
    return unique_sorted

def find_mrt(thresh, arr, time, a, n, state):
    times = []
    for i in range(n):
        if arr[i][a][40000] > thresh and state == 'upper':
            transit_time = check_thresh(arr[i][a], thresh, 0)
            times.append(time[transit_time])
        if arr[i][a][40000] < thresh and state == 'lower':
            transit_time = check_thresh(arr[i][a], thresh, 1)
            times.append(time[transit_time])        
    return times

def HS(x,x0,n,l):
    return (1+l*((x/x0)**n))/(1+((x/x0)**n))

def step(x, tau):
    return 1 * (x > tau)

def integration_const(p,time,time2,index,p_index,noise):
    dt = time[1] - time[0]
    points = time.size - 1

    st_10 = 0
    st_01 = 0
    st_00 = 0
    st_11 = 0
    
    D = np.empty(points + 1)
    E = np.empty(points + 1)

    D[0] = p['Dic'][index]
    E[0] = p['Eic'][index]

    lEtoD = p['lEtoD'][p_index]
    lDtoE = p['lDtoE'][p_index]
    
    gD = p['gD'][p_index]
    gE = p['gE'][p_index]
    
    kD = p['kD'][p_index]
    kE = p['kE'][p_index]
    
    trdEtoD = p['E0'][p_index]
    trdDtoE = p['D0'][p_index]
    
    nEtoD = p['nEtoD'][p_index]
    nDtoE = p['nDtoE'][p_index]
    
    tau = p['tau']
    
    thresh_D = (gD / kD) / 2 
    thresh_E = (gE / kE) / 2

    for i in range(1, points + 1):
        if time[i] % 1000 == 0:
            print(time[i])


        if time[i] in time2:
            lEtoD = max(0,min(1,(lEtoD + random.gauss(0,noise))))
            lDtoE = max(0,min(1,(lDtoE + random.gauss(0,noise))))
        

        D[i] = D[i-1] + dt * (gD * HS(E[int(i - 1 - p['tau'] * step(i, tau))], trdEtoD, nEtoD, lEtoD) - kD * D[i-1])
        E[i] = E[i-1] + dt * (gE * HS(D[int(i - 1 - p['tau'] * step(i, tau))], trdDtoE, nDtoE, lDtoE) - kE * E[i-1])
        
        
        if (time[i] > 50.0):
            if D[i] > thresh_D and E[i] > thresh_E:
                st_11 += 0.01
            elif D[i] > thresh_D and E[i] < thresh_E:
                st_10 += 0.01
            elif D[i] < thresh_D and E[i] > thresh_E:
                st_01 += 0.01
            elif D[i] < thresh_D and E[i] < thresh_E:
                st_00 += 0.01

    return D,E,st_10,st_01,st_11,st_00



T = 1000
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
dt2 = 0.01
time2 = np.arange(50.0,T+dt,dt2).round(2)
noises = np.arange(0.0, 0.0011,0.0002).round(4)
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
l_values = [[1,1]]#[0.1,0.1],[0.9,0.9],[0.1,0.9],[0.9,0.1]]
p['lAtoB'] = [0.025056, 0.057252, 0.012882, 0.019468, 0.021675]
p['lBtoC'] = [0.012458, 0.010081, 0.012339, 0.015197, 0.011492]
p['lCtoA'] = [0.010947, 0.012998, 0.050211, 0.029954, 0.013929]
p['lDtoE'] = [0.011075, 0.028405, 0.017531, 0.015665, 0.011402]
p['lEtoD'] = [0.021791, 0.018186, 0.018946, 0.020940, 0.018222]
#p['lDtoC'] = [0.1, 0.5, 0.9]
#p['lCtoD'] = [0.1, 0.5, 0.9]

#Hill Coeff
p['nAtoB'] = [2.000000, 6.000000, 3.000000, 2.000000, 3.000000]
p['nBtoC'] = [4.000000, 4.000000, 5.000000, 5.000000, 3.000000]
p['nCtoA'] = [6.000000, 4.000000, 4.000000, 5.000000, 5.000000]
p['nDtoE'] = [3.000000, 6.000000, 6.000000, 2.000000, 5.000000]
p['nEtoD'] = [5.000000, 3.000000, 3.000000, 5.000000, 6.000000]
p['nDtoC'] = 5.000000
p['nCtoD'] = 5.000000

taus = [0.0, 50.0, 100.0]

n = 25

for i in range(0, 1, 1):
    
    mrt_10 = []
    mrt_01 = []
    mrt_11 = []
    mrt_00 = []
    
    transit_times_lower = []

    if os.path.exists(folder + 'npy_files/size_50_random'+str(i)+'.npy'):
        arr = np.load(folder + 'npy_files/size_50_random'+str(i)+'.npy')
        p['Dic'] = arr[3]
        p['Eic'] = arr[4]
    else:
        #Initial Conditions
        p['Dic'] = np.random.random(size = n) * (1.5 * (p['gD'][i]/p['kD'][i]))
        p['Eic'] = np.random.random(size = n) * (1.5 * (p['gE'][i]/p['kE'][i]))

        arr = np.array([p['Aic'],p['Bic'],p['Cic'],p['Dic'],p['Eic']])

        np.save(folder + 'npy_files/size_50_random'+str(i)+'.npy', arr)
    
    for tau in range(len(taus)):
        mrts_upper = []
        stds_upper = []
        mrts_lower = []
        stds_lower = []
        p['tau'] = taus[tau]

        for noise in noises:
            arr_dynamics = []
            print(noise)
            
            states = {
                '10' : [],
                '01' : [],
                '11' : [],
                '00' : []
                }

            
            f, ax = plt.subplots(1,2,figsize = (10, 6), sharex = False, sharey = False) 
            
            for l in range(n):
                D,E,st_10,st_01,st_11,st_00 = integration_const(p, time, time2, l, i, noise)
                arr_dynamics.append(np.array([D,E]))
                
                states['10'].append(st_10)
                states['01'].append(st_01)
                states['00'].append(st_00)
                states['11'].append(st_11)
                
                ax[0].plot(time, D)
                ax[0].set_ylabel('D')
                
                ax[1].plot(time, E)
                ax[1].set_ylabel('E')
                
            f.tight_layout()
            plt.show()
            plt.close()
            
            mrt_10.append(mean(states['10']))
            mrt_01.append(mean(states['01']))
            mrt_11.append(mean(states['11']))
            mrt_00.append(mean(states['00']))
    
    #file = 'value_storage/toggleSwitch_transitTimesUpper_noise_dt2_' + str(dt2) + '_time_' + str(T) + '_param_' + str(i) + '.npy'
    #np.save(folder + file, np.array(transit_times_upper))
    
    #file = 'value_storage/toggleSwitch_transitTimesLower_noise_dt2_' + str(dt2) + '_time_' + str(T) + '_param_' + str(i) + '.npy'
    #np.save(folder + file, np.array(transit_times_lower))
    
   
    
f, ax = plt.subplots(2,2,figsize = (12, 12), sharex = False, sharey = False) 

for i in range(len(taus)):
    mrts_10 = []
    mrts_01 = []
    mrts_00 = []
    mrts_11 = []
    for j in range(len(noises)):
        index = j + i * len(noises)
        mrts_10.append(mrt_10[index])
        mrts_01.append(mrt_01[index])
        mrts_11.append(mrt_11[index])
        mrts_00.append(mrt_00[index])
        
    ax[0,0].plot(noises, mrts_10, label = str(taus[i]), marker = 'o')
    ax[0,0].set_ylabel('State 10')
    ax[0,0].legend()
    
    ax[0,1].plot(noises, mrts_01, label = str(taus[i]), marker = 'o')
    ax[0,1].set_ylabel('State 01')
    ax[0,1].legend()
    
    ax[1,0].plot(noises, mrts_00, label = str(taus[i]), marker = 'o')
    ax[1,0].set_ylabel('State 00')
    ax[1,0].legend()
    
    ax[1,1].plot(noises, mrts_11, label = str(taus[i]), marker = 'o')
    ax[1,1].set_ylabel('State 11')
    ax[1,1].legend()
        
f.tight_layout()
plt.xlabel('Noise')
plt.show()

            
                
            