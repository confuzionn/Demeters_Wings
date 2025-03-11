# ----------------------------------------------------------- Setup ----------------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette()

# ----------------------------------------------------------- Stall Speed ----------------------------------------------------------- #

# CL_max is within 1.3 - 1.9
xlim = 100
N=100
rho_cruise = 1.8685e-3 # slugs/ft^3, assuming average cruise of 8,000 ft
v_stall = 1.68781 * 100 # Convert kts to ft/s
CL_max = np.array([1.3, 1.5, 1.7, 1.9]) # Variable range of CL_max values based on Table 3.1 in Roskam

# Calculate wing loading based on stall speed
WS_stall = 0.5 * rho_cruise * v_stall**2 * CL_max # lb/ft^2 
print(WS_stall)

# Plot W/S vs W/P
plt.figure()
for i, ws in enumerate(WS_stall):
    plt.axvline(x=ws, color=colors[i], linestyle='--', alpha=0.7, label=f'$C_{{L_{{\\text{{max}}}}}}$ = {CL_max[i]:.1f}')
plt.title('W/P - W/S for Stall Speed')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("W/P $(lbf/hp)$")
plt.xlim(0, max(WS_stall)+5)
plt.ylim(0, 0.5)
plt.legend(loc='best')
plt.show()

# ----------------------------------------------------------- Takeoff ----------------------------------------------------------- #

S_TO_G = 1000 # ft
S_TO = 1500 # ft
rho_SL = 23.77e-4 #slugs/ft3
rho_takeoff = 20.48e-4 #slugs/ft3, assuming takeoff altitude of 5000 ft

# Relating S_TO to TOP_23
TOP_23 = (-8.134 + np.sqrt(8.314**2 - 4 * 0.0149 * -S_TO)) / (2 * 0.0149)
print('TOP23: {} '.format(TOP_23))

W_S = np.linspace(1,xlim,N) # Create W/S array
k_s = 1.2 # Takeoff speed factor
v_TO = k_s * v_stall # Takeoff velocity

# Plotting W/P vs W/S
plt.figure()
for CL_max in CL_max:
    WP_TO = (TOP_23 * (rho_takeoff / rho_SL) * CL_max) / (W_S)
    plt.plot(W_S, WP_TO, label=rf'$C_{{L_{{\text{{max}}, \text{{TO}}}}}} = {CL_max}$')
plt.title('W/P - W/S for Takeoff')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("W/P $(lbf/hp)$")
plt.xlim(0,60)
plt.ylim(0,100)
plt.legend(loc='best')
plt.show()

# ----------------------------------------------------------- Landing ----------------------------------------------------------- #

#Balanced field length
BFL= TOP_23*37.5
print(f"Required field length: {BFL} feet")

#Calculating landing distance
CL_max_values = np.array([1.3, 1.5, 1.7, 1.9])
W_TO = 19000 #lbs
W_L = 17000 #lbs
S_ref = 534.37 #lifting area
S_a = 600 # 600 ft for general aviation
for CL_max in CL_max_values:
    S_L_G = 80 * (W_L / S_ref) / (rho_takeoff/rho_SL * CL_max) + S_a
    print(f"For CL_max = {CL_max}: Landing ground roll distance = {S_L_G} feet")

#Calculate W/S
S_land = BFL *  0.6
W_S =  (rho_takeoff/rho_SL * CL_max_values) / (80*0.65) * (S_land - S_a)

# Plot W/S vs W/P
plt.figure()
for i, ws in enumerate(W_S):
    plt.axvline(x=ws, color=colors[i], linestyle='--', alpha=0.7, label=f'$C_{{L_{{\\text{{max}}}}}}$ = {CL_max_values[i]:.2f}')
plt.title('W/P - W/S for Landing')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("W/P $(lbf/hp)$")
plt.xlim(0, max(W_S)+5)
plt.ylim(0, 0.5)
plt.legend(loc='best')
plt.show()

#TW calculation
TW_landing = np.tile(W_S, 50)
print(f"TW_landing: {TW_landing}")

# ----------------------------------------------------------- Climb ----------------------------------------------------------- #

# FAR 23 Requirements

def PW_vs_WS(
        G,
        CD,
        CL,
        WS,
        rho
):
    PW = (G + (CL/CD)**(-1))/(CL**(1/2)) * (WS**(1/2))/(18.97*0.8*(rho/rho_SL))

    return PW

WS = np.linspace(1, xlim, N)

rho_5k = 20.48e-4 #slug/ft^3
CD0 = 0.01849
k = 0.04198
G_LVL1 = 0.083 # 8.3% for landplanes
ks = 1.2
# MTOW, flaps down, AEO, 0.94*T_TO

case1_PW = []

for i in WS:
    case1_PW.append((PW_vs_WS(G=G_LVL1,CD=CD0,CL=1.9,WS=i,rho=rho_SL)*((W_TO)/W_TO)**(3/2))**(-1))


G_LVL1_2 = 0.015 # at 5000ft
# MTOW, flaps up, OEI, T_TO

case2_PW = []

for i in WS:
    case2_PW.append((PW_vs_WS(G=G_LVL1_2,CD=CD0,CL=1.3,WS=i,rho=rho_5k)*(12/11)*((W_TO-2000)/W_TO)**(3/2))**(-1))



G_balk = 0.03
ks_balk = 1.3
#MTOW, landing flaps, AEO

case3_PW = []

for i in WS:
    case3_PW.append((PW_vs_WS(G=G_balk,CD=CD0,CL=1.9,WS=i,rho=rho_SL)*((W_TO-2000)/W_TO)**(3/2))**(-1))

plt.figure()
plt.plot(WS,case1_PW,label="Level 1, CL=1.9, AEO")
plt.plot(WS,case2_PW,label="Level 1, CL=1.3, OEI")
plt.plot(WS,case3_PW,label="Balked, CL=1.9, AEO")
plt.xlabel('W/S')
plt.ylabel('W/P')
plt.legend()
plt.title('W/P - W/S for Climb')
plt.show()

# ----------------------------------------------------------- Cruise ----------------------------------------------------------- #

def cruise_power_sizing(v_cruise, rho_cruise, eta_prop, CD0, k, xlim, N):
    """
    Compute the weight-to-power ratio (W/P) for cruise flight with a density-based power adjustment.

    Parameters:
    v_cruise: Cruise velocity in knots
    rho_cruise: Air density at cruise altitude (slugs/ft³)
    eta_prop: Propeller efficiency
    CD0: Zero-lift drag coefficient
    k: Induced drag constant (1 / (pi * e * AR))

    Returns:
    WP_cruise: Weight-to-power ratio (W/P) in lbf/hp
    WS: Wing loading (W/S) in lbf/ft²
    """
    
    rho_SL = 23.77e-4  # Sea level air density (slugs/ft³)
    WS = np.linspace(1, xlim, N)  # Wing loading values in lbf/ft²
    v_cruise = 1.6878 * v_cruise  # Convert knots to ft/s
    q_cruise = 0.5 * rho_cruise * v_cruise**2  # Dynamic pressure (psf)
    
    # Compute power loading W/P
    WP_cruise = 1 / (
        ((q_cruise * v_cruise * (CD0 + (k * WS**2) / (q_cruise**2))) / (550 * eta_prop * WS))
        * (rho_SL / rho_cruise) ** 0.75
    )
    
    return WP_cruise, WS

# Use function
WP_cruise, WS = cruise_power_sizing(v_cruise=150, 
                                rho_cruise=rho_cruise,
                                eta_prop=0.8, 
                                CD0=CD0, 
                                k=k,
                                xlim=xlim,
                                N=N
)

# Plot
plt.figure()
plt.plot(WS,WP_cruise,label="$v_{cruise}$ = 150 kts")
plt.xlabel('W/S $(lbf/ft^2)$')
plt.ylabel('W/P $(lbf/hp)$')
plt.title('W/P - W/S for Cruise @ 8000 ft')
plt.xlim(0,60)
plt.legend()
plt.show()

# ----------------------------------------------------------- Ceiling ----------------------------------------------------------- #

# Constants
W = 16856  # Takeoff Weight (lb)
S_ref = 654.375  # Reference Wing Area (ft^2)
e = 0.7558  # Oswald's Efficiency Factor
Cd0 = 0.01510  # Zero Lift Drag Coefficient
AR = 8.2556  # Aspect Ratio
prop_eff = 0.8  # Propeller Efficiency

max_alt = 25000  # Maximum Altitude (ft)
rho = 10.66e-4  # Air Density at 25,000 ft (lb/ft^3)
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
plt.title("W/P - W/S for Ceiling")
plt.legend()
plt.grid(True)
plt.show()

# ----------------------------------------------------------- Maneuver ----------------------------------------------------------- #

# Given parameters
CD_0 = 0.0151  # Zero-lift drag coefficient
n = 1 / np.cos(np.deg2rad(45))  # Load factor for 45° bank
eta_p= 0.8  # Propeller efficiency

# Define wing loading range
WS = np.linspace(0.1, xlim, N)

# Dynamic pressure (using V = 110 kts)
v_MN = 1.68781 * 110  # ft/s
q_MN = 0.5 * rho_SL * v_MN**2

R = 1 / (np.sqrt(n**2 - 1) * (32.17 / (v_MN**2)))
print(R)

TW_MN = (q_MN * CD_0) / WS + ((n/q_MN)**2 * WS * k)

# Compute W/P using velocity, prop efficiency, and HP conversion
WP_MN = (v_MN * eta_p) / (550 * TW_MN) # lbf/hp

# Plot
plt.figure(figsize=(8, 6))
plt.plot(WS, WP_MN, label="Maneuver Constraint", color='b')
plt.xlim(0, 60)
plt.xlabel("W/S $(lb/ft²)$")
plt.ylabel("W/P $(lb/hp)$")
plt.title("W/P - W/S for Maneuver")
plt.legend()
plt.grid()
plt.show()

# ----------------------------------------------------------- Total Constraint ----------------------------------------------------------- #
CL_max_19 = 1.9 # Setting CL_max = 1.9 for stall, takeoff, and landing

WS = np.linspace(1, xlim, N) # Create WS array
WS_stall = 0.5 * rho_cruise * v_stall**2 * CL_max_19
WP_TO_19 = (TOP_23 * (rho_takeoff / rho_SL) * CL_max_19) / (WS)
WS_landing =  (rho_takeoff/rho_SL * CL_max_19) / (80*0.65) * (S_land - S_a)

# Find intersection
intersection_idx = np.argmin(np.abs(WP_MN - WP_TO_19))  # Find index of minimum difference
intersection_ws = WS[intersection_idx]  # W/S value at intersection
intersection_wp = WP_MN[intersection_idx]  # W/P value at intersection
print('Design Parameters: {} $lbf/ft^2$, {} $lbf/hp$'.format(intersection_ws, intersection_wp))

plt.figure(figsize=(8,4),facecolor='#e3dfd7ff')
ax = plt.gca()  #
ax.set_facecolor("#e3dfd7ff")
plt.axvline(x=WS_stall, label='Stall', color='#284e3fff', linestyle='-', linewidth=2)
plt.plot(WS, WP_TO_19, label='Takeoff field length', linestyle='--', linewidth=2, color='#284e3fff')
plt.axvline(x=WS_landing, label='Landing field length', color='#284e3fff', linestyle='-', linewidth=2)
plt.plot(WS,case1_PW,label="Level 1 Climb, CL=1.9, AEO", linestyle='-.', color='#284e3fff')
plt.plot(WS,case2_PW,label="Level 1 Climb, CL=1.3, OEI", linestyle='-.',color='#284e3fff')
plt.plot(WS,case3_PW,label="Balked Climb, CL=1.9, AEO", linestyle='-.', color='#284e3fff')
plt.plot(WS,WP_cruise, label='Cruise', linestyle=':', linewidth=2, color='#284e3fff')
plt.plot(WS, WP_MN, label="Maneuver", linestyle='-', linewidth=2, color='#284e3fff')
plt.scatter(intersection_ws, intersection_wp, color='#cc0000ff', zorder=5, label='Design Point')
plt.fill_between(WS, 0, np.minimum(WP_MN, WP_TO_19), where=(WS <= WS_stall), color='#356854ff', alpha=0.5, label="Feasible Region")
plt.title('W/P - W/S', color='#301900ff')
plt.xlabel("W/S $(lb/ft^2)$", color='#301900ff')
plt.ylabel("W/P $(lb/hp)$", color='#301900ff')
plt.xlim(0,100)
plt.ylim(0,100)
plt.grid(True)
# plt.legend(loc='upper right')
legend = ax.legend(loc='upper right')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#e3dfd7ff')
plt.show()