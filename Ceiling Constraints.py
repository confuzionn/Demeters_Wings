import numpy as np
import matplotlib.pyplot as plt


S_ref = 534.375  # ft^2
e = 0.7558  # Oswald's Efficiency Factor
Cd0 = 0.01849  # Zero Lift Drag Coefficient
AR = 10.03268  # Aspect Ratio

max_alt = 25,000    #ft
rho = 0.0343    #lb/ft^3


k = 1 / (np.pi * e * AR)    #induced drag factor
T_W_min = 2 * np.sqrt(k * Cd0)

W_S_vals = np.linspace(0, 150, 100)
abs_ceiling_vals = [T_W_min] * 100

print(f"Induced Drag Factor: k = {k:.5f}")
print(f"Minimum Thrust over Weight: T/W = {T_W_min:.4f}")

plt.figure(figsize=(8, 6))
plt.plot(W_S_vals, abs_ceiling_vals, label="Absolute Ceiling", color='b')

plt.xlim(0, 150)
plt.ylim(0, T_W_min * 3)

plt.xlabel("W/S (lb / ft^2)")
plt.ylabel("T/W")
plt.title("Ceiling Design Constraint")
plt.legend()
plt.show()
