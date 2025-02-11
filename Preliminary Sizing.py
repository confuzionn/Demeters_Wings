import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# CL_max is within 1.3 - 1.9

rho = 0.002378 # kg/m^3
v_stall = 1.68781 * 100 # Convert kts to ft/s
CL_max = 1.3

W_S = 0.5 * rho * v_stall**2 * CL_max
print(W_S)

plt.figure()
plt.plot(W_S)
plt.axvline(x=W_S, color='red', linestyle='--')