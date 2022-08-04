import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
#from scipy import signal
from scipy.signal import find_peaks
from tqdm import tqdm
from itertools import groupby
from sklearn.metrics import auc

def runs(difference=0.011):
    start = None
    def inner(n):
        nonlocal start
        if start is None:
            start = n
        elif abs(start-n) > difference:
            start = n
        return start
    return inner

#print([next(g) for k, g in groupby(mylist, runs())])

folder = 'C:/project/Summer_2022_DrMKJolly/RS_TS/plots/dynamic/random_indi/delay/value_storage/'
file0 = 'delayDynamics_noise_0.0_dt2_0.5_time_10000_param_2_tauRange.npy'
file1 = 'delayDynamics_noise_0.0002_dt2_0.5_time_10000_param_2_tauRange.npy'
file2 = 'delayDynamics_noise_0.0004_dt2_0.5_time_10000_param_2_tauRange.npy'
file3 = 'delayDynamics_noise_0.0006_dt2_0.5_time_10000_param_2_tauRange.npy'
file4 = 'delayDynamics_noise_0.0008_dt2_0.5_time_10000_param_2_tauRange.npy'
file5 = 'delayDynamics_noise_0.001_dt2_0.5_time_10000_param_2_tauRange.npy'

dt=0.01
nic=10 ## number of initial conditions
zoom_factor=2
taus = np.arange(0.0, 101.0, 20.0)
noise = np.arange(0.0, 0.0012, 0.0002).round(4)
arr_peaks = {'peaks0': [], 'peaks1': [], 'peaks2': [], 'peaks3': [], 'peaks4': [], 'peaks5': []}
arr_auc = {'auc0': [], 'auc1': [], 'auc2': [], 'auc3': [], 'auc4': [], 'auc5': []}

for a,file in enumerate([file0, file1, file2, file3, file4, file5]):
    arr_dynamics_full = np.load(folder + file)
    for i in range(len(taus)):
        arr_dynamics = arr_dynamics_full[i]
        for l in tqdm(range(nic)):
            if l==0:
                N = arr_dynamics[l][0].size
                xf = fftfreq(N, dt)[50:N//zoom_factor];
                A_yf = 2.0/N * np.abs(fft(arr_dynamics[l][0])[50:N//zoom_factor])
                B_yf = 2.0/N * np.abs(fft(arr_dynamics[l][1])[50:N//zoom_factor])
                C_yf = 2.0/N * np.abs(fft(arr_dynamics[l][2])[50:N//zoom_factor])
            
            else: 
                A_yf += 2.0/N * np.abs(fft(arr_dynamics[l][0])[50:N//zoom_factor])
                B_yf += 2.0/N * np.abs(fft(arr_dynamics[l][1])[50:N//zoom_factor])
                C_yf += 2.0/N * np.abs(fft(arr_dynamics[l][2])[50:N//zoom_factor])
            
        A_yf=A_yf/nic; B_yf=B_yf/nic; C_yf=C_yf/nic;

        '''
        f, ax = plt.subplots(1,3,figsize = (16,8), sharex = True, sharey = True)
        ax[0].plot(xf, A_yf, linewidth = 1)
        ax[0].set_ylabel("A",fontsize=14)
        ax[0].grid()
        
        ax[1].plot(xf, B_yf, linewidth = 1)
        ax[1].set_ylabel("B",fontsize=14)
        ax[1].grid()
        
        ax[2].plot(xf, C_yf, linewidth = 1)
        ax[2].set_ylabel("C",fontsize=14)
        ax[2].grid()
        
        plt.tick_params(axis='y',labelsize=12,rotation=90)
        plt.tick_params(axis='x',labelsize=12)
        ax[1].set_xlabel("Freq (Hz)",fontsize=14)
        f.tight_layout()
        plt.show()
        '''
        
        print("Peak analysis for A")
        peaksA, _ = find_peaks(A_yf,prominence=0.5)
        uniq_peaksA = np.sort(np.unique(np.round(xf[peaksA],2)))
        uniq_peaksA = np.array([next(g) for k, g in groupby(uniq_peaksA, runs())])
        print (uniq_peaksA)
        #print(np.where(xf < 0.09 - 0.025)[0][-1])
        aucA = []
        
        for j in uniq_peaksA:
            x_range = xf[np.where(xf < j - 0.025)[0][-1] : np.where(xf > j + 0.025)[0][0]]
            y_range = A_yf[np.where(xf < j - 0.025)[0][-1] : np.where(xf > j + 0.025)[0][0]]
            aucA.append(auc(x_range, y_range).round(4) * 10000)
        print(aucA)
        
        print("Peak analysis for B")
        peaksB, _ = find_peaks(B_yf,prominence=0.5)
        uniq_peaksB = np.sort(np.unique(np.round(xf[peaksB],2)))
        uniq_peaksB = np.array([next(g) for k, g in groupby(uniq_peaksB, runs())])
        print (uniq_peaksB)        
        aucB = []
        
        for j in uniq_peaksB:
            x_range = xf[np.where(xf < j - 0.025)[0][-1] : np.where(xf > j + 0.025)[0][0]]
            y_range = B_yf[np.where(xf < j - 0.025)[0][-1] : np.where(xf > j + 0.025)[0][0]]
            aucB.append(auc(x_range, y_range).round(4) * 10000)
        print(aucB)
        
        print("Peak analysis for C")
        peaksC, _ = find_peaks(C_yf,prominence=0.5)
        uniq_peaksC = np.sort(np.unique(np.round(xf[peaksC],2)))
        uniq_peaksC = np.array([next(g) for k, g in groupby(uniq_peaksC, runs())])
        print (uniq_peaksC)
        aucC = []
        
        for j in uniq_peaksC:
            x_range = xf[np.where(xf < j - 0.025)[0][-1] : np.where(xf > j + 0.025)[0][0]]
            y_range = C_yf[np.where(xf < j - 0.025)[0][-1] : np.where(xf > j + 0.025)[0][0]]
            aucC.append(auc(x_range, y_range).round(4) * 10000)
        print(aucC)
        
        arr_peaks['peaks' + str(a)].append([uniq_peaksA, uniq_peaksB, uniq_peaksC])
        arr_auc['auc' + str(a)].append([aucA, aucB, aucC])
        
''' 
for a in range(6):
    f1, ax1 = plt.subplots(1,3,figsize = (16,8), sharex = True, sharey = True)
    for j in range(len(taus)):
        
        ax1[0].scatter(np.ones(len(arr_peaks['peaks' + str(a)][j][0])) * taus[j], arr_peaks['peaks' + str(a)][j][0], s = arr_auc['auc' + str(a)][j][0], color = 'red')
        ax1[0].set_ylabel("A",fontsize=14)
        
        ax1[1].scatter(np.ones(len(arr_peaks['peaks' + str(a)][j][1])) * taus[j], arr_peaks['peaks' + str(a)][j][1], s = arr_auc['auc' + str(a)][j][1], color = 'red')
        ax1[1].set_ylabel("B",fontsize=14)
        
        ax1[2].scatter(np.ones(len(arr_peaks['peaks' + str(a)][j][2])) * taus[j], arr_peaks['peaks' + str(a)][j][2], s = arr_auc['auc' + str(a)][j][2], color = 'red')
        ax1[2].set_ylabel("C",fontsize=14)
        
        
    
    ax1[0].grid()
    ax1[1].grid()
    ax1[2].grid()
    
    plt.tick_params(axis='y',labelsize=12,rotation=90)
    plt.tick_params(axis='x',labelsize=12)
    ax1[1].set_xlabel("Taus",fontsize=14)
    f1.tight_layout()
'''

for j in range(len(taus)):
    f1, ax1 = plt.subplots(1,3,figsize = (24 ,8), sharex = True, sharey = True)
    for k in range(len(noise)):
        
        ax1[0].scatter((np.ones(len(arr_peaks['peaks' + str(k)][j][0])) * noise[k]).astype('str'), arr_peaks['peaks' + str(k)][j][0], s = arr_auc['auc' + str(k)][j][0], color = 'red')
        ax1[0].set_ylabel("A",fontsize=14)
        
        ax1[1].scatter((np.ones(len(arr_peaks['peaks' + str(k)][j][1])) * noise[k]).astype('str'), arr_peaks['peaks' + str(k)][j][1], s = arr_auc['auc' + str(k)][j][1], color = 'red')
        ax1[1].set_ylabel("B",fontsize=14)
        
        ax1[2].scatter((np.ones(len(arr_peaks['peaks' + str(k)][j][2])) * noise[k]).astype('str'), arr_peaks['peaks' + str(k)][j][2], s = arr_auc['auc' + str(k)][j][2], color = 'red')
        ax1[2].set_ylabel("C",fontsize=14)
        
        
    
    ax1[0].grid()
    ax1[1].grid()
    ax1[2].grid()
    
    plt.tick_params(axis='y',labelsize=12,rotation=90)
    plt.tick_params(axis='x',labelsize=12)
    ax1[1].set_xlabel("Noise",fontsize=14)
    plt.title(str(taus[j]), fontsize = 14)
    f1.tight_layout()
