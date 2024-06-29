import numpy as np
import matplotlib.pyplot as plt

def x(n, l, H0, g, k):
    ans = n / ( l * H0 / ( g / k ))#x = n/(λ*H0/(g/k))
    return ans

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


for i in range(4,6,1):
    xAtoB = x(p['nAtoB'][i], p['lAtoB'][i], p['trdAtoB'][i], p['gA'][i], p['kA'][i])
    xAtoC = x(p['nAtoC'][i], p['lAtoC'][i], p['trdAtoC'][i], p['gA'][i], p['kA'][i])
    xBtoA = x(p['nBtoA'][i], p['lBtoA'][i], p['trdBtoA'][i], p['gB'][i], p['kB'][i])
    xBtoC = x(p['nBtoC'][i], p['lBtoC'][i], p['trdBtoC'][i], p['gB'][i], p['kB'][i])
    xCtoA = x(p['nCtoA'][i], p['lCtoA'][i], p['trdCtoA'][i], p['gC'][i], p['kC'][i])
    xCtoB = x(p['nCtoB'][i], p['lCtoB'][i], p['trdCtoB'][i], p['gC'][i], p['kC'][i])
    
    x_axis = ['A--|B / B--|A', 'C--|A / A--|C', 'B--|C / C--|B', 'B--|A / A--|B', 'A--|C / C--|A', 'C--|B / B--|C']
    y_axis_unnorm = [xAtoB/xBtoA, xCtoA/xAtoC, xBtoC/xCtoB, xBtoA/xAtoB, xAtoC/xCtoA, xCtoB/xBtoC]
    y_axis = [float(i)/sum(y_axis_unnorm) for i in y_axis_unnorm]
    
    plt.bar(x_axis, y_axis)
    
    plt.tick_params(axis='x', labelrotation=10, labelsize= 10)
    plt.xlabel('Link',fontsize=15)
    plt.ylabel('Link Strength',fontsize=15)
    plt.title('Param'+str(i),fontsize=15)
    plt.tight_layout()
    plt.show()
    plt.close()

