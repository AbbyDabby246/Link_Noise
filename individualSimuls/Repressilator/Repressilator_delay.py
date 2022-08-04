import numpy as np
import os
import random

def HS(x,x0,n,l):
    return (1+l*((x/x0)**n))/(1+((x/x0)**n))

def step(x, tau):
    return 1 * (x > tau)

def integration_const(p,time,time2,index,p_index,noise):
    dt = time[1] - time[0]
    points = time.size - 1
    #dt2 = time2[1] - time2[0]

    #lamda_random = []

    A = np.empty(points + 1)
    B = np.empty(points + 1)
    C = np.empty(points + 1)
    #D = np.empty(points + 1)
    #E = np.empty(points + 1)

    A[0] = p['Aic'][index]
    B[0] = p['Bic'][index]
    C[0] = p['Cic'][index]
    #D[0] = p['Dic'][index]
    #E[0] = p['Eic'][index]

    #lEtoD = p['lEtoD'][p_index]
    #lDtoE = p['lDtoE'][p_index]
    lAtoB = p['lAtoB'][p_index]
    lBtoC = p['lBtoC'][p_index]
    lCtoA = p['lCtoA'][p_index]

    for i in range(1, points + 1):
        if time[i] % 1000 == 0:
            print(time[i])

        # Introducing Noise in the fold change values
        if time[i] in time2:
            #lEtoD = max(0,min(1,(lEtoD + random.gauss(0,noise))))#max(lDtoC + random.uniform(0.1,0.5),0) if random.random()<0.5 else max(lDtoC - random.uniform(0.1,0.5),0)
            #lDtoE = max(0,min(1,(lDtoE + random.gauss(0,noise))))
            lAtoB = max(0,min(1,(lAtoB + random.gauss(0,noise))))
            lBtoC = max(0,min(1,(lBtoC + random.gauss(0,noise))))
            lCtoA = max(0,min(1,(lCtoA + random.gauss(0,noise))))#max(lCtoD + random.uniform(0.1,0.5),0) if random.random()<0.5 else max(lCtoD - random.uniform(0.1,0.5),0)
            #print(time[i], lEtoD, lDtoE,lAtoB,lBtoC,lCtoA)
        #lamda_random.append([lDtoE,lDtoE,lAtoB,lBtoC,lCtoA])

        #solving using eulers method
        A[i] = A[i-1] + dt * (p['gA'][p_index] * HS(C[int(i - 1 - p['tau'] * step(i, p['tau']))], p['C0'][p_index], p['nCtoA'][p_index], lCtoA) - p['kA'][p_index] * A[i-1])
        B[i] = B[i-1] + dt * (p['gB'][p_index] * HS(A[int(i - 1 - p['tau'] * step(i, p['tau']))], p['A0'][p_index], p['nAtoB'][p_index], lAtoB) - p['kB'][p_index] * B[i-1])
        C[i] = C[i-1] + dt * (p['gB'][p_index] * HS(B[int(i - 1 - p['tau'] * step(i, p['tau']))], p['B0'][p_index], p['nBtoC'][p_index], lBtoC) - p['kC'][p_index] * C[i-1])

        #D[i] = D[i-1] + dt * (p['gD'][p_index] * HS(E[int(i - 1 - p['tau'] * step(i, p['tau']))], p['E0'][p_index], p['nEtoD'][p_index], lEtoD) - p['kD'][p_index] * D[i-1])
        #E[i] = E[i-1] + dt * (p['gE'][p_index] * HS(D[int(i - 1 - p['tau'] * step(i, p['tau']))], p['D0'][p_index], p['nDtoE'][p_index], lDtoE) - p['kE'][p_index] * E[i-1])

    return A,B,C#,D,E,lamda_random



T = 10000
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
dt2 = 0.5
time2 = np.arange(1000.0,9000+dt,dt2).round(2)
noises = [0.0006] #variance of normal dist
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
l_values = [[1,1]]
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

taus = np.arange(0.0, 101.0, 20.0)

n = 10

#combs = [[l_values[i],l_values[j]] for i in range(4) for j in range(4) if i!=j]


for i in range(2, 3, 1):

    if os.path.exists(folder + 'npy_files/size_50_random'+str(i)+'.npy'):
        arr = np.load(folder + 'npy_files/size_50_random'+str(i)+'.npy')
        p['Aic'] = arr[0]
        p['Bic'] = arr[1]
        p['Cic'] = arr[2]
        #p['Dic'] = arr[3]
        #p['Eic'] = arr[4]
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
        arr_dynamics = []
        print(noise)
        
        for tau in taus:
            p['tau'] = tau
            tau_dynamics = []
            
            for l in range(n):
                A,B,C = integration_const(p, time, time2, l, i, noise)
                tau_dynamics.append(np.array([A,B,C]))
            
            arr_dynamics.append(tau_dynamics)
                
        file = 'value_storage/delayDynamics_noise_' + str(noise) + '_dt2_' + str(dt2) + '_time_' + str(T) + '_param_' + str(i) + '_tauRange.npy'
        np.save(folder + file, np.array(arr_dynamics))
