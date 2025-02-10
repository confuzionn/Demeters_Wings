import numpy as np

S_ref = 534.375     #ft^2
e = 0.7558      #Oswald's Efficiency Factor
CdO = 0.01849   #Zero Lift Drag Coefficient
AR = 10.03268   #Aspect Ratio

k = 1 / (np.pi * e * AR)

T_W_ratio_min = 2 * np.sqrt(k * Cd0)