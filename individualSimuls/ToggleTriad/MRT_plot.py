import matplotlib.pyplot as plt
import numpy as np
import re

noises = np.arange(0.0, 0.00051, 0.0001).round(4)
states = {
          '100' : [],
          '010' : [],
          '001' : [],
          '110' : [],
          '011' : [],
          '101' : [],
          '111' : []
            }
param = 3
file = 'stateMRTs_param' + str(param) + '.txt'
f = open('C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleTriad/MRTs/' + file, 'r')

for row in f:
  row = re.split(r'\s{1,}', row)
  row = row[:-1]
  if str(row[0]) in list(states.keys()):
    print(row)
    states[row[0]].append(round(float(row[-1]),2))
f.close()


f = plt.figure()
f.set_figwidth(7)
f.set_figheight(5)
for key in states.keys():
    y = states[key]
    x = noises
    plt.plot(x,y, label = key)
plt.title('MRTs for TT; Param: ' + str(param))
plt.yscale('symlog')
plt.legend()
plt.show()