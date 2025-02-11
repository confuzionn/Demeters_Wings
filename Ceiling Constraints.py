import numpy as np
import matplotlib.pyplot as plt

# Constants
W = 16856  # Takeoff Weight (lb)
S_ref = 534.375  # Reference Wing Area (ft^2)
e = 0.7558  # Oswald's Efficiency Factor
Cd0 = 0.01849  # Zero Lift Drag Coefficient
AR = 10.03268  # Aspect Ratio
prop_eff = 0.8  # Propeller Efficiency
CL_clean_1 = 1.3  # Lift Coefficient (Clean Configuration 1)
CL_clean_2 = 1.9  # Lift Coefficient (Clean Configuration 2)

max_alt = 25000  # Maximum Altitude (ft)
rho = 0.0343  # Air Density at 25,000 ft (lb/ft^3)

# Induced Drag Factor
k = 1 / (np.pi * e * AR)

# Minimum Thrust-to-Weight Ratio
T_W_min = 2 * np.sqrt(k * Cd0)

# Wing Loading Values (start from 0.1 to avoid division by zero)
W_S_vals = np.linspace(0.1, 150, 100)

# Absolute Ceiling Values
abs_ceiling_vals = [T_W_min] * 100

# Velocity and W/P Calculations for CL = 1.3
V_1 = np.sqrt((2 * W) / (rho * CL_clean_1))
W_P_1 = (T_W_min ** (-1)) * prop_eff / V_1
W_P_1_vals = [W_P_1] * 100

# Velocity and W/P Calculations for CL = 1.9
V_2 = np.sqrt((2 * W) / (rho * CL_clean_2))
W_P_2 = (T_W_min ** (-1)) * prop_eff / V_2
W_P_2_vals = [W_P_2] * 100
# Print Results
print(f"Induced Drag Factor: k = {k:.5f}")
print(f"Minimum Thrust over Weight: T/W = {T_W_min:.4f}")

# Plot 1: Absolute Ceiling (T/W vs W/S)
plt.figure(figsize=(8, 6))
plt.plot(W_S_vals, abs_ceiling_vals, label="Absolute Ceiling", color='b')

plt.xlim(0, 150)
plt.ylim(0, T_W_min * 3)

plt.xlabel("W/S (lb / ft^2)")
plt.ylabel("T/W")
plt.title("Ceiling Design Constraint")
plt.legend()
plt.grid(True)
plt.show()

# Plot 2: Weight-to-Power Ratio (W/P vs W/S)
plt.figure(figsize=(8, 6))
plt.plot(W_S_vals, W_P_1_vals, label="W/P (CL = 1.3)", color='r', linestyle='--', alpha=0.7)
plt.plot(W_S_vals, W_P_2_vals, label="W/P (CL = 1.9)", color='g', linestyle='-.', alpha=0.7)

plt.xlim(0, 150)
plt.ylim(0, 0.1)  # Adjust y-axis limit dynamically

plt.xlabel("W/S (lb / ft^2)")
plt.ylabel("W/P (lb/hp)")
plt.title("Weight-to-Power Ratio Constraint")
plt.legend()
plt.grid(True)
plt.show()
