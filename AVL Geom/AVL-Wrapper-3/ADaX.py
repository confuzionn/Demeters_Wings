from Code.Codebase import *
from Aircraft_Geometry import Update_Geometry
from AVL_Presets import *
import numpy as np
import matplotlib.pyplot as plt

#Update_Geometry()
ad = ADaX('AR5_stick.avl')
# ad.input('mass avl/planes/bd.mass')
# ad.input('case avl/planes/bd.run')
# ad.input('mode')
# ad.input('n')
# ad.input('w cameron_test.eig')
ad.run()
Clx = 1.25

Stall_Prediction(ad, Clx)

case = 1

eigen_data = np.loadtxt('AR5_stick.eig',
            skiprows=(case-1)*8+3,
            max_rows=8,
            usecols=[1, 2]
           )

nat_freq = np.sqrt(np.square(eigen_data[:,0])+np.square(eigen_data[:,1]))/(2*np.pi)
damp_rat = -eigen_data[:,0]/(nat_freq*2*np.pi)
period = 2*np.pi/(eigen_data[:,1])
t_half = 0.693/(eigen_data[:,0])

save_freqy = []
save_freqx = []
save_nat = []
save_laty = []
save_latx = []

for i in range(0,8):
    if np.abs(eigen_data[i,1]) > 0:
        save_freqy.append(eigen_data[i,1])
        save_freqx.append(eigen_data[i,0])
        save_nat.append(nat_freq[i])
    else :
        save_laty.append(eigen_data[i,1])
        save_latx.append(eigen_data[i,0])

ind_freq = np.argsort(save_nat)

print("Natural Frequency:", nat_freq)
print("Damping Ratio:", damp_rat)
print("Period:", period)
print("T 1/2:", t_half)

yoffset = 0.2
xoffset = 1
plt.scatter(eigen_data[:,0],eigen_data[:,1])

plt.text(min(save_latx)-xoffset,min(save_laty)+yoffset,'Roll')

plt.text(max(save_latx)-xoffset,max(save_laty)+yoffset,'Spiral')

plt.text(save_freqx[ind_freq[-1]]-xoffset,save_freqy[ind_freq[-1]]+yoffset,'Short')
plt.text(save_freqx[ind_freq[-1-1]]-xoffset,save_freqy[ind_freq[-1-1]]+yoffset,'Short')

plt.text(save_freqx[ind_freq[2]]-xoffset,save_freqy[ind_freq[2]]+yoffset,'Dutch Roll')
plt.text(save_freqx[ind_freq[3]]-xoffset,save_freqy[ind_freq[3]]+yoffset,'Dutch Roll')

plt.text(save_freqx[ind_freq[0]]-xoffset,save_freqy[ind_freq[0]]+yoffset,'Phugoid')
plt.text(save_freqx[ind_freq[1]]-xoffset,save_freqy[ind_freq[1]]+yoffset,'Phugoid')



plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.title('Eigenvalues for Case {0}'.format(case))
plt.grid()
plt.show()