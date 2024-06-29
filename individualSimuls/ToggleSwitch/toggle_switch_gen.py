import numpy as np
import os
import random
import matplotlib.pyplot as plt
from statistics import stdev, mean
import pandas as pd
#import multiprocessing as mp

np.random.seed(123456)
no_of_cpu_for_use=5

def check_thresh(l,thr,temp = 0):
    if temp:
        c = [n > thr for n in l[50000:]]
    else:
        c = [n < thr for n in l[50000:]]
    
    if sum(c) > 0:
        return c.index(1)
    return -1

def count_switches(lst):
    switch_count = 0
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1]:
            switch_count += 1
    return switch_count

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

def integration_snapshots(p,time,time2,index,p_index,noise, nSnaps = 50):
    dt = time[1] - time[0]
    points = time.size - 1

    st_10 = 0
    st_01 = 0
    st_00 = 0
    st_11 = 0

    curr_st = []
    
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

    lamda_track = [[lEtoD],[lDtoE]]
    for i in range(1, points + 1):
        if time[i] % 1000 == 0:
            print(time[i])
          

        if time[i] == 100:
            lEtoD = max(0,min(1,(lEtoD + random.gauss(0,noise))))
            lDtoE = max(0,min(1,(lDtoE + random.gauss(0,noise))))
        
        lamda_track[0].append(lEtoD)
        lamda_track[1].append(lDtoE) 

        D[i] = D[i-1] + dt * (gD * HS(E[int(i - 1 - p['tau'] * step(i, tau))], trdEtoD, nEtoD, lEtoD) - kD * D[i-1])
        E[i] = E[i-1] + dt * (gE * HS(D[int(i - 1 - p['tau'] * step(i, tau))], trdDtoE, nDtoE, lDtoE) - kE * E[i-1])
        
        
        if (time[i] > 50.0):
            if D[i] > thresh_D and E[i] > thresh_E:
                st_11 += 0.01
                curr_st.append('st_11')
            elif D[i] > thresh_D and E[i] < thresh_E:
                st_10 += 0.01
                curr_st.append('st_10')
            elif D[i] < thresh_D and E[i] > thresh_E:
                st_01 += 0.01
                curr_st.append('st_01')
            elif D[i] < thresh_D and E[i] < thresh_E:
                st_00 += 0.01
                curr_st.append('st_00')

    num_switches = count_switches(curr_st)

    return D,E,st_10,st_01,st_11,st_00, lamda_track, num_switches


def integration_const(p,time,time2,index,p_index,noise):
    dt = time[1] - time[0]
    points = time.size - 1

    st_10 = 0
    st_01 = 0
    st_00 = 0
    st_11 = 0

    curr_st = []
    
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

    lamda_track = [[lEtoD],[lDtoE]]
    for i in range(1, points + 1):
        if time[i] % 1000 == 0:
            print(time[i])
          

        if time[i] in time2:
            lEtoD = max(0,min(1,(lEtoD + random.gauss(0,noise))))
            lDtoE = max(0,min(1,(lDtoE + random.gauss(0,noise))))
        
        lamda_track[0].append(lEtoD)
        lamda_track[1].append(lDtoE) 

        D[i] = D[i-1] + dt * (gD * HS(E[int(i - 1 - p['tau'] * step(i, tau))], trdEtoD, nEtoD, lEtoD) - kD * D[i-1])
        E[i] = E[i-1] + dt * (gE * HS(D[int(i - 1 - p['tau'] * step(i, tau))], trdDtoE, nDtoE, lDtoE) - kE * E[i-1])
        
        
        if (time[i] > 50.0):
            if D[i] > thresh_D and E[i] > thresh_E:
                st_11 += 0.01
                curr_st.append('st_11')
            elif D[i] > thresh_D and E[i] < thresh_E:
                st_10 += 0.01
                curr_st.append('st_10')
            elif D[i] < thresh_D and E[i] > thresh_E:
                st_01 += 0.01
                curr_st.append('st_01')
            elif D[i] < thresh_D and E[i] < thresh_E:
                st_00 += 0.01
                curr_st.append('st_00')

    num_switches = count_switches(curr_st)

    return D,E,st_10,st_01,st_11,st_00, lamda_track, num_switches



T = 1000
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
dt2 = 0.01
time2 = np.arange(50.0,T+dt,dt2).round(2)
noises = np.arange(0.0, 0.11,0.02).round(4)
n_points = time.size
print(time2)

cont_noise = pd.DataFrame(columns = ['Param', 'Noise', 'Run', 'D', 'E'])
noncont_noise = pd.DataFrame(columns = ['Param', 'Noise', 'Run', 'D', 'E'])

folder = "C:/project/Summer_2022_DrMKJolly/RS_TS/plots/dynamic/random_indi/delay/"

p = {}

#production paramters
p['gD'] = [91.947511, 50.212433, 64.102184, 83.975461, 39.691404]
p['gE'] = [38.696894, 95.256533, 25.889225, 66.912087, 77.538958]

#Threshold constants
p['D0'] = [53.920926, 14.888436, 31.018099, 53.370193,  3.534087]
p['E0'] = [36.738347, 35.491987,  7.507220, 10.120514, 42.641427]

#Degradation constants
p['kD'] = [0.699470, 0.390561, 0.829410, 0.207175, 0.881824]
p['kE'] = [0.442699, 0.444709, 0.801409, 0.498639, 0.431325]

#Fold Change
p['lDtoE'] = [0.011075, 0.028405, 0.017531, 0.015665, 0.011402]
p['lEtoD'] = [0.021791, 0.018186, 0.018946, 0.020940, 0.018222]

#Hill Coeff
p['nDtoE'] = [3.000000, 6.000000, 6.000000, 2.000000, 5.000000]
p['nEtoD'] = [5.000000, 3.000000, 3.000000, 5.000000, 6.000000]

taus = [0.0]

n = 25

for i in range(0, 5, 1):
    lamda_track = dict()
    
    mrt_10 = []
    mrt_01 = []
    mrt_11 = []
    mrt_00 = []
    
    transit_times_lower = []

    switching_events = {}

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
            
            lamda_track[noise] = []

            states = {
                '10' : [],
                '01' : [],
                '11' : [],
                '00' : []
                }
            
            total_switches = 0

            
            #f, ax = plt.subplots(1,2,figsize = (10, 6), sharex = False, sharey = False) 
            f = plt.figure()
            f.set_figwidth(7)
            f.set_figheight(5)
            for l in range(n):
                D,E,st_10,st_01,st_11,st_00,lt, num_sw = integration_const(p, time, time2, l, i, noise)
                D_snap,E_snap,st_10_snap,st_01_snap,st_11_snap,st_00_snap, lt_snap, num_sw_snap = integration_snapshots(p, time, time2, l, i, noise)
                arr_dynamics.append(np.array([D,E]))

                cont_noise = pd.concat([cont_noise, pd.DataFrame({'Param': i, 'Noise': noise, 'Run': l, 'D': D, 'E': E})], ignore_index = True, axis = 1)
                noncont_noise = pd.concat([noncont_noise, pd.DataFrame({'Param': i, 'Noise': noise, 'Run': l, 'D': D_snap, 'E': E_snap})], ignore_index = True)

                total_switches += num_sw

                plt.plot(D_snap, E_snap)
                plt.title('nonCont A vs B; Param: ' + str(i) + ' Noise: ' + str(noise))
                 

                lamda_track[noise].append(lt_snap)
                
                states['10'].append([st_10, st_10_snap])
                states['01'].append([st_01, st_01_snap])
                states['00'].append([st_00, st_00_snap])
                states['11'].append([st_11, st_11_snap])
                
                #ax[0].plot(time, D)
                #ax[0].set_ylabel('D')
                
                #ax[1].plot(time, E)
                #ax[1].set_ylabel('E')
            
            #f.tight_layout()
            #plt.show()
            #plt.close()
            plt.savefig("snap_AvsB_Island_Noise_" + str(noise) + "_param_" + str(i) +'.png')
            plt.close()

            switching_events[noise] = total_switches

            mrt_10.append([mean([i[0] for i in states['10']]), mean([i[1] for i in states['10']])])
            mrt_01.append([mean([i[0] for i in states['01']]), mean([i[1] for i in states['01']])])
            mrt_11.append([mean([i[0] for i in states['11']]), mean([i[1] for i in states['11']])])
            mrt_00.append([mean([i[0] for i in states['00']]), mean([i[1] for i in states['00']])])
    
    f = plt.figure()
    f.set_figwidth(7)
    f.set_figheight(5)
    
    plt.bar(range(len(switching_events)), list(switching_events.values()), align='center')
    plt.xticks(range(len(switching_events)), list(switching_events.keys()))
    plt.title("Num switching evenets; NOISE: " + str(noise)+ ", Param: " + str(i))

    #plt.savefig("Num_switching_NOISE_" + str(noise)+ "_param_" + str(i) +'.png')
    plt.close()

    f, ax = plt.subplots(ncols=2, sharex = False, sharey = False) 
    #f = plt.figure()
    f.set_figwidth(15)
    f.set_figheight(5)

    for k in range(len(taus)):
        mrts_10 = []
        mrts_01 = []
        mrts_00 = []
        mrts_11 = []

        mrts_10_snap = []
        mrts_01_snap = []
        mrts_00_snap = []
        mrts_11_snap = []
        for j in range(len(noises)):
            index = j + k * len(noises)
            mrts_10.append(mrt_10[index][0])
            mrts_01.append(mrt_01[index][0])
            mrts_11.append(mrt_11[index][0])
            mrts_00.append(mrt_00[index][0])

            mrts_10_snap.append(mrt_10[index][1])
            mrts_01_snap.append(mrt_01[index][1])
            mrts_00_snap.append(mrt_00[index][1])
            mrts_11_snap.append(mrt_11[index][1])
            
        #ax[0,0].plot(noises, mrts_10, label = '10', marker = 'o')
        #ax[0,0].set_ylabel('State 10')
        #ax[0,0].legend()
        
        #ax[0,1].plot(noises, mrts_01, label = '01', marker = 'o')
        #ax[0,1].set_ylabel('State 01')
        #ax[0,1].legend()
        
        #ax[1,0].plot(noises, mrts_00, label = '00', marker = 'o')
        #ax[1,0].set_ylabel('State 00')
        #ax[1,0].legend()
        
        #ax[1,1].plot(noises, mrts_11, label = '11', marker = 'o')
        #ax[1,1].set_ylabel('State 11')
        #ax[1,1].legend()
        
        ax[0].plot(noises, mrts_10, label='10',marker='o', alpha = 0.7)
        ax[0].plot(noises, mrts_01, label='01',marker='o', alpha = 0.7)
        ax[0].plot(noises, mrts_11, label='11',marker='o', alpha = 0.7)
        ax[0].plot(noises, mrts_00, label='00',marker='o', alpha = 0.7)

        ax[0].set_title("Cont Noise")
        ax[0].set_yscale("symlog")
        ax[0].set_ylabel('MRT')
        ax[0].set_xlabel('Noise')
        ax[0].legend()

        ax[1].plot(noises, mrts_10_snap, label='10',marker='o', alpha = 0.7)
        ax[1].plot(noises, mrts_01_snap, label='01',marker='o', alpha = 0.7)
        ax[1].plot(noises, mrts_11_snap, label='11',marker='o', alpha = 0.7)
        ax[1].plot(noises, mrts_00_snap, label='00',marker='o', alpha = 0.7)

        ax[1].set_title("Snapshots")
        ax[1].set_yscale("symlog")
        ax[1].set_ylabel('MRT')
        ax[1].set_xlabel('Noise')
        ax[1].legend()
            
    #f.tight_layout()
    #plt.xlabel('Noise')
    #plt.show()
    #plt.close()
    f.suptitle('Symlog MRTs for TS; Param: ' + str(i))
    #plt.yscale('symlog')
    #plt.legend()
    plt.savefig('combined_MRTs_sym_Param_' + str(i) + '.png')
    plt.close()

    
    f = plt.figure()
    f.set_figwidth(7)
    f.set_figheight(5)

    for noise in noises:
        for num in range(5):
            plt.plot(time, lamda_track[noise][num][0])
            plt.plot(time, lamda_track[noise][num][1])

    plt.title('Snapshot lamda track; Param: ' + str(i))
    plt.savefig('snap_lamda_track_Param_' + str(i) + '.png')
    #plt.yscale('symlog')
    #plt.legend()
    plt.close()  

#cont_noise.to_csv('ToggleSwitch_contDynamics.csv')
#noncont_noise.to_csv('ToggleSwitch_noncontDynamics.csv')
    
print(cont_noise.head())
print("\n")
print(noncont_noise.head())



            
                
            