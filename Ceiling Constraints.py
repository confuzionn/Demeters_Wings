import numpy as np
import matplotlib.pyplot as plt

# Constants
W = 16856  # Takeoff Weight (lb)
S_ref = 654.375  # Reference Wing Area (ft^2)
e = 0.7558  # Oswald's Efficiency Factor
Cd0 = 0.01510  # Zero Lift Drag Coefficient
AR = 8.2556  # Aspect Ratio
prop_eff = 0.8  # Propeller Efficiency

max_alt = 25000  # Maximum Altitude (ft)
rho = 0.0343  # Air Density at 25,000 ft (lb/ft^3)
R_C_max = 100 #ft/min

# Induced Drag Factor
k = 1 / (np.pi * e * AR)

CL = np.sqrt(Cd0 / k)

# Minimum Thrust-to-Weight Ratio
T_W_min = 2 * np.sqrt(k * Cd0)
T_W_serv = T_W_min + R_C_max *((Cd0 / k) ** (1/4)) * ((rho / 2) ** (1/2)) * ((W / S_ref) ** (-1/2))

# Wing Loading Values (start from 0.1 to avoid division by zero)
W_S_vals = np.linspace(0.1, 150, 100)

# Absolute Ceiling Values
abs_ceiling_vals = [T_W_min] * 100

# Velocity and W/P Calculations for CL = 1.3
V = 1.68781 * 110  # ft/s

W_P = V * prop_eff / (550 * T_W_min)
W_P_vals = [W_P] * 100
# Print Results
print(f"Induced Drag Factor: k = {k:.5f}")
print(f"Minimum Thrust over Weight: T/W = {T_W_min:.4f}")

plt.figure(figsize=(8, 6))
plt.plot(W_S_vals, W_P_vals, label="W/P (CL = 0.66)", color='r', linestyle='--', alpha=0.7)

plt.xlim(0, 150)
plt.ylim(0, 10)  # Adjust y-axis limit dynamically

plt.xlabel("W/S (lb / ft^2)")
plt.ylabel("W/P (lb/hp)")
plt.title("Weight-to-Power Ratio Constraint")
plt.legend()
plt.grid(True)
plt.show()
