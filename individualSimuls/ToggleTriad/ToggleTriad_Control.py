# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:30:14 2022

@author: abhay
"""

import numpy as np
import os
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import statistics as stat

def z_scores(tt):
    A_mean = stat.mean([j[5900] for i in tt for j in i])
    A_std  = stat.stdev([j[5900] for i in tt for j in i])
    A_z = [round((i[0][5900] - A_mean) / A_std, 2) for i in tt]
    
    B_mean = stat.mean([i[1][5900] for i in tt])
    B_std  = stat.stdev([i[1][5900] for i in tt])
    B_z = [round((i[1][5900] - A_mean) / A_std, 2) for i in tt]
    
    C_mean = stat.mean([i[2][5900] for i in tt])
    C_std  = stat.stdev([i[2][5900] for i in tt])
    C_z = [round((i[2][5900] - A_mean) / A_std, 2) for i in tt]
    
    return [A_z, B_z, C_z]


def H(b,b0a,nba):
    return 1 / (1 + (b / b0a) ** nba)

def HS(b, b0a, nba, lba):
    return H(b, b0a, nba) + lba * (1 - H(b, b0a, nba))

def step(x, tau):
    return 1 * (x > tau)

def integration_const(p,time,index,p_index):
    dt = time[1] - time[0]
    points = time.size - 1

    A = np.empty(points + 1)
    B = np.empty(points + 1)
    C = np.empty(points + 1)
    
    A[0] = p['Aic'][index]
    B[0] = p['Bic'][index]
    C[0] = p['Cic'][index]

    lAtoB = p['lAtoB'][p_index]
    lBtoC = p['lBtoC'][p_index]
    lCtoA = p['lCtoA'][p_index]
    lAtoC = p['lAtoC'][p_index]
    lCtoB = p['lCtoB'][p_index]
    lBtoA = p['lBtoA'][p_index]
    
    ga = p['gA'][p_index]
    gb = p['gB'][p_index]
    gc = p['gC'][p_index]
    
    ka = p['kA'][p_index]
    kb = p['kB'][p_index]
    kc = p['kC'][p_index]
    
    trdCtoA = p['trdCtoA'][p_index]
    trdBtoA = p['trdBtoA'][p_index]
    trdAtoB = p['trdAtoB'][p_index]
    trdCtoB = p['trdCtoB'][p_index]
    trdAtoC = p['trdAtoC'][p_index]
    trdBtoC = p['trdBtoC'][p_index]
    
    nCtoA = p['nCtoA'][p_index]
    nBtoA = p['nBtoA'][p_index]
    nAtoB = p['nAtoB'][p_index]
    nCtoB = p['nCtoB'][p_index]
    nAtoC = p['nAtoC'][p_index]
    nBtoC = p['nBtoC'][p_index]

    for i in range(1, points + 1):
        if time[i] % 100 == 0:
            print(time[i])

        '''
        if time[i] in time2:
            lAtoB = max(0,min(1,(lAtoB + random.gauss(0,noise))))
            lBtoC = max(0,min(1,(lBtoC + random.gauss(0,noise))))
            lCtoA = max(0,min(1,(lCtoA + random.gauss(0,noise))))
        '''

        
        A[i] = A[i-1] + dt * (ga * HS(C[i - 1], trdCtoA, nCtoA, lCtoA) * HS(B[i - 1], trdBtoA, nBtoA, lBtoA) - ka * A[i-1])
        B[i] = B[i-1] + dt * (gb * HS(A[i - 1], trdAtoB, nAtoB, lAtoB) * HS(C[i - 1], trdCtoB, nCtoB, lCtoB) - kb * B[i-1])
        C[i] = C[i-1] + dt * (gc * HS(B[i - 1], trdBtoC, nBtoC, lBtoC) * HS(A[i - 1], trdAtoC, nAtoC, lAtoC) - kc * C[i-1])

    return A,B,C



T = 100
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
n_points = time.size

folder = "C:/project/Summer_2022_DrMKJolly/TT/"

p = {}

#production paramters
p['gA'] = [12.018249, 5.230785, 92.301214, 59.528363, 43.125467, 74.661273]
p['gB'] = [2.675641, 76.596182, 47.940154, 79.558099, 42.705912, 88.839696]
p['gC'] = [60.869333, 6.714682, 46.216528, 89.527126, 22.316457, 33.342705]

#Threshold constants
p['trdAtoB'] = [8.091984 , 3.091889 , 0.400201, 2.367581, 7.893864, 5.979105 ]
p['trdBtoC'] = [2.561299 , 7.347394 ,  8.68131, 5.630365, 9.875833, 11.509482]
p['trdCtoA'] = [0.329027 , 6.569456 , 5.929507, 7.94269 , 9.151694, 4.745698 ]
p['trdAtoC'] = [9.705534 , 9.334428 , 5.253287, 4.583083, 6.194162, 9.969809 ]
p['trdCtoB'] = [4.111351 , 6.053961 , 7.564845, 3.941795, 12.18131, 8.414342 ]
p['trdBtoA'] = [11.313256, 10.165773, 6.924967, 7.630631, 9.723269, 11.243532]

#Degradation constants
p['kA'] = [0.804511, 0.577791, 0.361671, 0.940077, 0.631683, 0.7806  ]
p['kB'] = [0.178982, 0.534397, 0.949621, 0.509631, 0.758171, 0.913511]
p['kC'] = [0.767847, 0.180412, 0.230333, 0.899618, 0.429697, 0.678669]

#Fold Change
p['lAtoB'] = [0.080769, 0.013618, 0.117883, 0.012118, 0.013122, 0.020298]
p['lBtoC'] = [0.01015 , 0.030059, 0.019566, 0.014407, 0.013806, 0.022222]
p['lCtoA'] = [0.013822, 0.010811, 0.014476, 0.02025 , 0.017554, 0.012852]
p['lAtoC'] = [0.013274, 0.021265, 0.012013, 0.010072, 0.011343, 0.010311]
p['lCtoB'] = [0.017927, 0.031853, 0.011715, 0.011788, 0.069397, 0.016738]
p['lBtoA'] = [0.017287, 0.05123 , 0.022806, 0.012885, 0.024531, 0.056011]

#Hill Coeff
p['nAtoB'] = [6.000000, 4.000000, 1.000000, 6.000000, 4.000000, 6.000000]
p['nBtoC'] = [4.000000, 2.000000, 4.000000, 3.000000, 4.000000, 4.000000]
p['nCtoA'] = [1.000000, 5.000000, 4.000000, 6.000000, 5.000000, 5.000000]
p['nAtoC'] = [4.000000, 1.000000, 6.000000, 3.000000, 3.000000, 5.000000]
p['nCtoB'] = [5.000000, 4.000000, 3.000000, 2.000000, 6.000000, 2.000000]
p['nBtoA'] = [4.000000, 3.000000, 4.000000, 3.000000, 6.000000, 6.000000]

n = 100

for i in range(3, 4, 1):

    if os.path.exists(folder + 'value_storage/initial_conditions/size_' + str(n) + '_random' + str(i) + '.npy'):
        arr = np.load(folder + 'value_storage/initial_conditions/size_' + str(n) + '_random' + str(i) + '.npy')
        p['Aic'] = arr[0]
        p['Bic'] = arr[1]
        p['Cic'] = arr[2]
    else:
        #Initial Conditions
        p['Aic'] = np.random.random(size = n) * (1.5 * (p['gA'][i]/p['kA'][i]))
        p['Bic'] = np.random.random(size = n) * (1.5 * (p['gB'][i]/p['kB'][i]))
        p['Cic'] = np.random.random(size = n) * (1.5 * (p['gC'][i]/p['kC'][i]))

        arr = np.array([p['Aic'],p['Bic'],p['Cic']])

        np.save(folder + 'value_storage/initial_conditions/size_' + str(n) + '_random' + str(i) + '.npy', arr)

    
    tt_dynamics = []
    
    f, ax = plt.subplots(1,4,figsize = (20, 6), sharex = True, sharey = False)
    
    for l in range(n):
        A,B,C = integration_const(p, time, l, i)
        
        ax[0].plot(time[500:], A[500:], linewidth = 1.5, label = str(p['Aic'][l]))
        ax[0].set_ylabel("A",fontsize=14)
        ax[1].plot(time[500:], B[500:], linewidth = 1.5, label = str(p['Bic'][l]))
        ax[1].set_ylabel("B",fontsize=14)
        ax[2].plot(time[500:], C[500:], linewidth = 1.5, label = str(p['Cic'][l]))
        ax[2].set_ylabel("C",fontsize=14)
        
        ax[3].plot(time[500:], A[500:], linewidth = 1.5, label = 'A', color = 'red', alpha = 0.1)
        ax[3].plot(time[500:], B[500:], linewidth = 1.5, label = 'B', color = 'blue', alpha = 0.1)
        ax[3].plot(time[500:], C[500:], linewidth = 1.5, label = 'C', color = 'green', alpha = 0.1)
        ax[3].set_ylabel("ABC",fontsize=14)
        
        tt_dynamics.append(np.array([A,B,C]))
        
    plt.tick_params(axis='y',labelsize=12,rotation=90)
    plt.tick_params(axis='x',labelsize=12)
    f.suptitle('TT_Control')
    ax[1].set_xlabel("Time (in weeks)",fontsize=14)
    
    red_patch = mpatches.Patch(color='red', label='A')
    blue_patch = mpatches.Patch(color='blue', label='B')
    green_patch = mpatches.Patch(color='green', label='C')
    
    ax[3].legend(handles = [red_patch, blue_patch, green_patch])
    
    f.tight_layout()
    plt.show()
    #plt.savefig('plot_'+'lDtoC_lCtoD_' + str(lamdas[0]) + '_' +str(lamdas[1])+'_random_'+str(dt2)+'.jpeg',dpi=1000)
    plt.close()
    
    #arr_dynamics.append(tt_dynamics)
        
    file = 'value_storage/toggleTriad_control_time_' + str(T) + '_param_' + str(i) + '.npy'
    np.save(folder + file, np.array(tt_dynamics))
    
    z_score = z_scores(tt_dynamics)
    
    print('Z_scores for A: ', list(set(z_score[0])))
    print('Z_scores for B: ', list(set(z_score[1])))
    print('Z_scores for C: ', list(set(z_score[2])))
    
    f, ax = plt.subplots(1,3,figsize = (15, 6), sharex = True, sharey = False)
    
    ax[0].hist(z_score[0], bins = 4)
    ax[0].set_ylabel("A_freq",fontsize=14)
    ax[1].hist(z_score[1], bins = 4)
    ax[1].set_ylabel("B_freq",fontsize=14)
    ax[2].hist(z_score[2], bins = 4)
    ax[2].set_ylabel("C_freq",fontsize=14)
    
    f.tight_layout()
    plt.show()
    plt.close()