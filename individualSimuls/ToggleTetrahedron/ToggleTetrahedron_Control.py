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

    D_mean = stat.mean([i[2][5900] for i in tt])
    D_std  = stat.stdev([i[2][5900] for i in tt])
    D_z = [round((i[2][5900] - A_mean) / A_std, 2) for i in tt]
    
    return [A_z, B_z, C_z, D_z]


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
    D = np.empty(points + 1)
    
    A[0] = p['Aic'][index]
    B[0] = p['Bic'][index]
    C[0] = p['Cic'][index]
    D[0] = p['Dic'][index]

    lAtoB = p['lAtoB'][p_index]
    lBtoC = p['lBtoC'][p_index]
    lCtoA = p['lCtoA'][p_index]
    lAtoC = p['lAtoC'][p_index]
    lCtoB = p['lCtoB'][p_index]
    lBtoA = p['lBtoA'][p_index]
    lAtoD = p['lAtoD'][p_index]
    lBtoD = p['lBtoD'][p_index]
    lCtoD = p['lCtoD'][p_index]
    lDtoC = p['lDtoC'][p_index]
    lDtoB = p['lDtoB'][p_index]
    lDtoA = p['lDtoA'][p_index]
    
    ga = p['gA'][p_index]
    gb = p['gB'][p_index]
    gc = p['gC'][p_index]
    gd = p['gD'][p_index]
    
    ka = p['kA'][p_index]
    kb = p['kB'][p_index]
    kc = p['kC'][p_index]
    kd = p['kD'][p_index]
    
    trdCtoA = p['trdCtoA'][p_index]
    trdBtoA = p['trdBtoA'][p_index]
    trdAtoB = p['trdAtoB'][p_index]
    trdCtoB = p['trdCtoB'][p_index]
    trdAtoC = p['trdAtoC'][p_index]
    trdBtoC = p['trdBtoC'][p_index]
    trdCtoD = p['trdCtoD'][p_index]
    trdBtoD = p['trdBtoD'][p_index]
    trdAtoD = p['trdAtoD'][p_index]
    trdDtoB = p['trdDtoB'][p_index]
    trdDtoC = p['trdDtoC'][p_index]
    trdDtoA = p['trdDtoA'][p_index]
    
    nCtoA = p['nCtoA'][p_index]
    nBtoA = p['nBtoA'][p_index]
    nAtoB = p['nAtoB'][p_index]
    nCtoB = p['nCtoB'][p_index]
    nAtoC = p['nAtoC'][p_index]
    nBtoC = p['nBtoC'][p_index]
    nCtoD = p['nCtoD'][p_index]
    nBtoD = p['nBtoD'][p_index]
    nAtoD = p['nAtoD'][p_index]
    nDtoB = p['nDtoB'][p_index]
    nDtoC = p['nDtoC'][p_index]
    nDtoA = p['nDtoA'][p_index]

    for i in range(1, points + 1):
        if time[i] % 100 == 0:
            print(time[i])

        '''
        if time[i] in time2:
            lAtoB = max(0,min(1,(lAtoB + random.gauss(0,noise))))
            lBtoC = max(0,min(1,(lBtoC + random.gauss(0,noise))))
            lCtoA = max(0,min(1,(lCtoA + random.gauss(0,noise))))
        '''

        
        A[i] = A[i-1] + dt * (ga * HS(C[i - 1], trdCtoA, nCtoA, lCtoA) * HS(B[i - 1], trdBtoA, nBtoA, lBtoA) * HS(D[i - 1], trdDtoA, nDtoA, lDtoA) - ka * A[i-1])
        B[i] = B[i-1] + dt * (gb * HS(A[i - 1], trdAtoB, nAtoB, lAtoB) * HS(C[i - 1], trdCtoB, nCtoB, lCtoB) * HS(D[i - 1], trdDtoB, nDtoB, lDtoB) - kb * B[i-1])
        C[i] = C[i-1] + dt * (gc * HS(B[i - 1], trdBtoC, nBtoC, lBtoC) * HS(A[i - 1], trdAtoC, nAtoC, lAtoC) * HS(D[i - 1], trdDtoC, nDtoC, lDtoC) - kc * C[i-1])
        D[i] = D[i-1] + dt * (gd * HS(A[i - 1], trdAtoD, nAtoD, lAtoD) * HS(B[i - 1], trdBtoD, nBtoD, lBtoD) * HS(C[i - 1], trdCtoD, nCtoD, lCtoD) - kd * D[i-1])

    return A,B,C,D



T = 100
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
n_points = time.size

folder = "C:/project/Summer_2022_DrMKJolly/TTetra/"

p = {}



#production paramters
p['gA'] = [49.250974, 21.291053, 11.855185]
p['gB'] = [19.394058, 33.94156 , 96.093647]
p['gC'] = [87.801864, 37.822804, 14.395325]
p['gD'] = [51.248971, 51.241041, 2.468148 ]

#Threshold constants
p['trdBtoA'] = [0.967474, 0.283163, 0.908393]
p['trdCtoA'] = [0.240827, 0.640190, 0.096357]
p['trdDtoA'] = [0.594364, 0.698889, 1.084239]
p['trdAtoB'] = [0.844277, 1.048565, 1.111341]
p['trdCtoB'] = [1.165158, 0.823650, 0.279576]
p['trdDtoB'] = [0.688072, 0.592589, 1.200326]
p['trdAtoC'] = [0.291655, 1.018860, 1.252788]
p['trdBtoC'] = [0.09921 , 0.020724, 1.080284]
p['trdDtoC'] = [0.295686, 0.145776, 0.788254]
p['trdAtoD'] = [0.99736 , 1.134829, 0.720232]
p['trdBtoD'] = [0.590067, 0.821745, 0.513065]
p['trdCtoD'] = [1.247305, 1.125203, 0.551829]

#Degradation constants
p['kA'] = [0.798911, 0.127533, 0.130997]
p['kB'] = [0.16937 , 0.587455, 0.95451 ]
p['kC'] = [0.592715, 0.720125, 0.6069 ]
p['kD'] = [0.146741, 0.571025, 0.610782]

#Fold Change
p['lBtoA'] = [0.012752, 0.017283, 0.351338]
p['lCtoA'] = [0.022607, 0.039878, 0.027123]
p['lDtoA'] = [0.050241, 0.010254, 0.016484]
p['lAtoB'] = [0.010739, 0.011151, 0.016373]
p['lCtoB'] = [0.015581, 0.015275, 0.010395]
p['lDtoB'] = [0.01099 , 0.021597, 0.012045]
p['lAtoC'] = [0.011552, 0.012367, 0.01354 ]
p['lBtoC'] = [0.167756, 0.010857, 0.044065]
p['lDtoC'] = [0.474966, 0.013967, 0.021406]
p['lAtoD'] = [0.195663, 0.01324 , 0.015594]
p['lBtoD'] = [0.056791, 0.014933, 0.024771]
p['lCtoD'] = [0.022647, 0.01647 , 0.013834]

#Hill Coeff
p['nBtoA'] = [6.000000, 7, 6 ]
p['nCtoA'] = [6.000000, 8, 10]
p['nDtoA'] = [10.00000, 6, 10]
p['nAtoB'] = [10.00000, 6, 8 ]
p['nCtoB'] = [8.000000, 9, 10]
p['nDtoB'] = [10.00000, 7, 7 ]
p['nAtoC'] = [9.000000, 9, 9 ]
p['nBtoC'] = [9.000000, 7, 6 ]
p['nDtoC'] = [7.000000, 7, 6 ]
p['nAtoD'] = [8.000000, 7, 8 ]
p['nBtoD'] = [8.000000, 9, 10]
p['nCtoD'] = [10.00000, 10, 8]

n = 5

for i in range(2, 3, 1):

    if os.path.exists(folder + 'value_storage/initial_conditions/size_' + str(n) + '_random' + str(i) + '.npy'):
        arr = np.load(folder + 'value_storage/initial_conditions/size_' + str(n) + '_random' + str(i) + '.npy')
        p['Aic'] = arr[0]
        p['Bic'] = arr[1]
        p['Cic'] = arr[2]
        p['Dic'] = arr[3]
    else:
        #Initial Conditions
        p['Aic'] = np.random.random(size = n) * (1.5 * (p['gA'][i]/p['kA'][i]))
        p['Bic'] = np.random.random(size = n) * (1.5 * (p['gB'][i]/p['kB'][i]))
        p['Cic'] = np.random.random(size = n) * (1.5 * (p['gC'][i]/p['kC'][i]))
        p['Dic'] = np.random.random(size = n) * (1.5 * (p['gD'][i]/p['kD'][i]))

        arr = np.array([p['Aic'],p['Bic'],p['Cic'],p['Dic']])

        np.save(folder + 'value_storage/initial_conditions/size_' + str(n) + '_random' + str(i) + '.npy', arr)

    
    tt_dynamics = []
    
    f, ax = plt.subplots(1,5,figsize = (25, 6), sharex = True, sharey = False)
    
    for l in range(n):
        A,B,C,D = integration_const(p, time, l, i)
        
        ax[0].plot(time[500:], A[500:], linewidth = 1.5, label = str(p['Aic'][l]))
        ax[0].set_ylabel("A",fontsize=14)
        ax[1].plot(time[500:], B[500:], linewidth = 1.5, label = str(p['Bic'][l]))
        ax[1].set_ylabel("B",fontsize=14)
        ax[2].plot(time[500:], C[500:], linewidth = 1.5, label = str(p['Cic'][l]))
        ax[2].set_ylabel("C",fontsize=14)        
        ax[3].plot(time[500:], D[500:], linewidth = 1.5, label = str(p['Dic'][l]))
        ax[3].set_ylabel("D",fontsize=14)
        
        ax[4].plot(time[500:], A[500:], linewidth = 1.5, label = 'A', color = 'red', alpha = 0.1)
        ax[4].plot(time[500:], B[500:], linewidth = 1.5, label = 'B', color = 'blue', alpha = 0.1)
        ax[4].plot(time[500:], C[500:], linewidth = 1.5, label = 'C', color = 'green', alpha = 0.1)
        ax[4].plot(time[500:], D[500:], linewidth = 1.5, label = 'D', color = 'purple', alpha = 0.1)
        ax[4].set_ylabel("ABCD",fontsize=14)
        
        tt_dynamics.append(np.array([A,B,C,D]))
        
    plt.tick_params(axis='y',labelsize=12,rotation=90)
    plt.tick_params(axis='x',labelsize=12)
    f.suptitle('TT_Control')
    ax[1].set_xlabel("Time (in weeks)",fontsize=14)
    
    red_patch = mpatches.Patch(color='red', label='A')
    blue_patch = mpatches.Patch(color='blue', label='B')
    green_patch = mpatches.Patch(color='green', label='C')    
    yellow_patch = mpatches.Patch(color='purple', label='D')
    
    ax[4].legend(handles = [red_patch, blue_patch, green_patch, yellow_patch])
    
    f.tight_layout()
    plt.show()
    #plt.savefig('plot_'+'lDtoC_lCtoD_' + str(lamdas[0]) + '_' +str(lamdas[1])+'_random_'+str(dt2)+'.jpeg',dpi=1000)
    plt.close()
    
    #arr_dynamics.append(tt_dynamics)

    """    
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
    """