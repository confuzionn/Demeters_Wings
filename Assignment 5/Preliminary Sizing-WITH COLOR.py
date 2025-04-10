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

'''
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
'''
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
'''
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
'''
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
'''
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
'''
#TW calculation
TW_landing = np.tile(W_S, 50)

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
'''
plt.figure()
plt.plot(WS,case1_PW,label="Level 1, CL=1.9, AEO")
plt.plot(WS,case2_PW,label="Level 1, CL=1.3, OEI")
plt.plot(WS,case3_PW,label="Balked, CL=1.9, AEO")
plt.xlabel('W/S')
plt.ylabel('W/P')
plt.legend()
plt.title('W/P - W/S for Climb')
'''
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
'''
# Plot
plt.figure()
plt.plot(WS,WP_cruise,label="$v_{cruise}$ = 150 kts")
plt.xlabel('W/S $(lbf/ft^2)$')
plt.ylabel('W/P $(lbf/hp)$')
plt.title('W/P - W/S for Cruise @ 8000 ft')
plt.xlim(0,60)
plt.legend()
'''
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

'''
plt.figure(figsize=(8, 6))
plt.plot(W_S_vals, W_P_vals, label="W/P (CL = 0.66)", color='r', linestyle='--', alpha=0.7)

plt.xlim(0, 150)
plt.ylim(0, 10)  # Adjust y-axis limit dynamically

plt.xlabel("W/S (lb / ft^2)")
plt.ylabel("W/P (lb/hp)")
plt.title("W/P - W/S for Ceiling")
plt.legend()
plt.grid(True)
'''
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

TW_MN = (q_MN * CD_0) / WS + ((n/q_MN)**2 * WS * k)

# Compute W/P using velocity, prop efficiency, and HP conversion
WP_MN = (v_MN * eta_p) / (550 * TW_MN) # lbf/hp

'''
# Plot
plt.figure(figsize=(8, 6))
plt.plot(WS, WP_MN, label="Maneuver Constraint", color='b')
plt.xlim(0, 60)
plt.xlabel("W/S $(lb/ft²)$")
plt.ylabel("W/P $(lb/hp)$")
plt.title("W/P - W/S for Maneuver")
plt.legend()
plt.grid()
'''

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
plt.axvline(x=WS_stall, label='Stall', color='k',  linestyle='-', linewidth=2)
plt.plot(WS, WP_TO_19, label='Takeoff field length', linestyle='--', linewidth=2)
plt.axvline(x=WS_landing, label='Landing field length', color='#284e3fff', linestyle='-', linewidth=2)
plt.plot(WS,case1_PW,label="Level 1 Climb, CL=1.9, AEO", linestyle='-.')
plt.plot(WS,case2_PW,label="Level 1 Climb, CL=1.3, OEI", linestyle='-.')
plt.plot(WS,case3_PW,label="Balked Climb, CL=1.9, AEO", linestyle='-.')
plt.plot(WS,WP_cruise, label='Cruise', linestyle=':', linewidth=2)
plt.plot(WS, WP_MN, label="Maneuver", linestyle='-', linewidth=2)
plt.scatter(intersection_ws, intersection_wp, color='red', zorder=5, label='Design Point')
plt.fill_between(WS, 0, np.minimum(WP_MN, WP_TO_19), where=(WS <= WS_stall), color='#356854ff', alpha=0.5, label="Feasible Region")
plt.title('W/P - W/S', color='#301900ff')
plt.xlabel("W/S $(lb/ft^2)$", color='#301900ff')
plt.ylabel("W/P $(lb/hp)$", color='#301900ff')
plt.xlim(0,100)
plt.ylim(0,100)
plt.grid(True)
legend = ax.legend(loc='upper right')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#e3dfd7ff')


# ----------------------------------------------------------- POWER BALANCE EQUATIONS ----------------------------------------------------------- #

from sympy import symbols, Eq, solve

eff_GT = 0.4
eff_PM = 0.99
eff_EM = 0.96
eff_P = 0.85
eff_GB = 0.96

def get_Powers(shaft_ratio, supplied_ratio, P_p):
    # Define symbols for the unknowns
    P_gt, P_f, P_s1, P_gb, P_p1, P_e1, P_e2, P_s2, P_p2, P_bat = symbols(
        'P_gt P_f P_s1 P_gb P_p1 P_e1 P_e2 P_s2 P_p2 P_bat')

    # Define the equations
    eq1 = Eq(P_gt, eff_GT * P_f)
    eq2 = Eq(P_s1 + P_gb, eff_GB * P_gt)
    eq3 = Eq(P_p1, eff_P * P_s1)
    eq4 = Eq(P_e1, eff_EM * P_gb)
    eq5 = Eq(P_e2, eff_PM * (P_e1 + P_bat))
    eq6 = Eq(P_s2, eff_EM * P_e2)
    eq7 = Eq(P_p2, eff_P * P_s2)
    eq8 = Eq(shaft_ratio, P_s2 / (P_s1 + P_s2))
    eq9 = Eq(supplied_ratio, P_bat / (P_bat + P_f))
    eq10 = Eq(P_p, P_p2 + P_p1)

    # Solve the system of equations
    solution = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10],
                     (P_gt, P_f, P_s1, P_gb, P_p1, P_e1, P_e2, P_s2, P_p2, P_bat))

    # Evaluate numerical values
    numerical_solution = {var: value.evalf() for var, value in solution.items()}

    return numerical_solution


# ----------------------------------------------------------- Electric Battery Power Loading ----------------------------------------------------------- #

WP_TO_bat = []
WP_MN_bat = []

for P_p in 1 / WP_TO_19:
    sol = get_Powers(1, 1, P_p)
    WP_TO_bat.append(1 / sol[symbols('P_bat')])

for P_p in 1 / WP_MN:
    sol = get_Powers(1, 1, P_p)
    WP_MN_bat.append(1 / sol[symbols('P_bat')])

# Ensure WP_MN_bat and WP_TO_bat are NumPy arrays
WP_TO_bat = np.array([float(val.evalf()) for val in WP_TO_bat])
WP_MN_bat = np.array([float(val.evalf()) for val in WP_MN_bat])

# Find index of minimum difference
intersection_idx2 = np.argmin(np.abs(WP_MN_bat - WP_TO_bat))

# W/S and W/P values at intersection
intersection_ws2 = WS[intersection_idx2]  # Assuming WS is a NumPy array or list with the same length as WP_MN_bat
intersection_wp2 = np.abs(WP_MN_bat[intersection_idx2])  # W/P value at intersection

# Print the result
print('Max W/P Design Parameters: {} lbf/ft^2, {} lbf/hp'.format(intersection_ws2, intersection_wp2))

WS = np.array(WS)

# Electric Battery Power Loading
plt.figure(figsize=(8,4),facecolor='#e3dfd7ff')
ax = plt.gca()  #
ax.set_facecolor("#e3dfd7ff")
plt.axvline(x=WS_stall, label='Stall', color='k', linestyle='-', linewidth=2)
plt.plot(WS, WP_TO_bat, label='Takeoff field length', linestyle='-', linewidth=2)
plt.plot(WS, WP_MN_bat, label="Maneuver", linestyle='-', linewidth=2)
plt.scatter(intersection_ws2, intersection_wp2, color='red', zorder=5, label='Design for min bat weight')
plt.fill_between(
    WS,
    0,
    np.minimum(WP_TO_bat, WP_MN_bat),
    where=(WS <= WS_stall),
    color='#284e3fff',
    alpha=0.5,
    label="Feasible Region"
)
plt.title('$W/P_{bat} - W/S$ for Electric Configuration')
plt.xlabel("$W/S (lb/ft^2)$",color='#301900ff')
plt.ylabel("$W/P_{bat} (lb/hp)$",color='#301900ff')
plt.xlim(0,100)
plt.ylim(0,100)
plt.grid(True)
legend = ax.legend(loc='upper right')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#e3dfd7ff')

W_elec = 12224 # lbf (from weight estimation code)
S_ref_elec = W_elec / WS
P_elec = W_elec / intersection_wp2
S_elec = W_elec / intersection_ws2
P_MN_elec = W_elec / WP_TO_bat
P_TO_elec = W_elec / WP_MN_bat
P_case1_elec =  W_elec/np.array(case1_PW)
P_case2_elec =  W_elec/np.array(case2_PW)
P_case3_elec = W_elec/np.array(case3_PW)
P_cruise_elec = W_elec/np.array(WP_cruise)
S_stall = W_elec / WS_stall
print('Electric P & S Design Parameters: {} ft^2, {} hp'.format(S_elec, P_elec))

# Electric Power vs Wing Area Dimensional Constraint
plt.figure(figsize=(8,4),facecolor='#e3dfd7ff')
ax = plt.gca()  #
ax.set_facecolor("#e3dfd7ff")
plt.axvline(x=S_stall, label='Stall', color='k', linestyle='-', linewidth=2)
plt.plot(S_ref_elec, np.abs(P_TO_elec), label='Takeoff field length', linestyle='-', linewidth=2)
plt.plot(S_ref_elec, P_case1_elec,label="Level 1 Climb, CL=1.9, AEO")
plt.plot(S_ref_elec, P_case2_elec,label="Level 1 Climb, CL=1.3, OEI")
plt.plot(S_ref_elec, P_case3_elec,label="Balked Climb, CL=1.9, AEO")
plt.plot(S_ref_elec, P_cruise_elec, label='Cruise', linestyle='-', linewidth=2)
plt.plot(S_ref_elec, np.abs(P_MN_elec), label="Maneuver", linestyle='-', linewidth=2)
plt.scatter(495, 1550, color='red', zorder=5, label='Design for min battery size')
plt.fill_between(S_ref_elec, np.maximum(P_TO_elec, P_MN_elec), 4000, where=(S_ref_elec >= S_stall), color='#284e3fff', alpha=0.5, label="Feasible Region")
plt.title('$P - S$ for Electric Configuration')
plt.xlabel("$S (ft^2)$")
plt.ylabel("$P (hp)$")
plt.xlim(200,600)
plt.ylim(0,4000)
plt.grid(True)
legend = ax.legend(loc='upper right')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#e3dfd7ff')

# ----------------------------------------------------------- Hybrid Battery Power Loading ----------------------------------------------------------- #

WP_TO_bat2 = []
WP_MN_bat2 = []

for P_p in 1 / WP_TO_19:
    sol = get_Powers(0.5, 0.5, P_p)
    WP_TO_bat2.append(1 / sol[symbols('P_bat')])

for P_p in 1 / WP_MN:
    sol = get_Powers(0.5, 0.5, P_p)
    WP_MN_bat2.append(1 / sol[symbols('P_bat')])

# Ensure WP_MN_gt and WP_TO_gt are NumPy arrays
WP_TO_bat2 = np.array([float(val.evalf()) for val in WP_TO_bat2])
WP_MN_bat2 = np.array([float(val.evalf()) for val in WP_MN_bat2])

# Find index of minimum difference
intersection_idx3 = np.argmin(np.abs(WP_MN_bat2 - WP_TO_bat2))

# W/S and W/P values at intersection
intersection_ws3 = WS[intersection_idx3]  # W/S value at intersection
intersection_wp3 = np.abs(WP_MN_bat2[intersection_idx3])  # W/P value at intersection

# Print the result
print('Max W/P Design Parameters: {} lbf/ft^2, {} lbf/hp'.format(intersection_ws3, intersection_wp3))

# Hybrid Power Loading
plt.figure(figsize=(8,4),facecolor='#e3dfd7ff')
ax = plt.gca()  #
ax.set_facecolor("#e3dfd7ff")
plt.axvline(x=WS_stall, label='Stall', color='k', linestyle='-', linewidth=2)
plt.plot(WS, WP_TO_bat2, label='Takeoff field length', linestyle='-', linewidth=2)
plt.plot(WS, WP_MN_bat2, label="Maneuver", linestyle='-', linewidth=2)
plt.scatter(intersection_ws3, intersection_wp3, color='red', zorder=5, label='Design point')
plt.fill_between(
    WS,
    0,
    np.minimum(WP_TO_bat2, WP_MN_bat2),
    where=(WS <= WS_stall),
    color='#284e3fff',
    alpha=0.5,
    label="Feasible Region"
)
plt.title('$W/P_{bat} - W/S$ for Hybrid Configuration')
plt.xlabel("$W/S (lb/ft^2)$")
plt.ylabel("$W/P_{bat} (lb/hp)$")
plt.xlim(0,100)
plt.ylim(0,100)
plt.grid(True)
legend = ax.legend(loc='upper right')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#e3dfd7ff')

W_hybrid = 7700 # lbf (from weight estimation code)
S_ref_hybrid = W_hybrid / WS
P_hybrid = W_hybrid / intersection_wp3
S_hybrid = W_hybrid / intersection_ws3
P_MN_hybrid = W_hybrid / WP_TO_bat2
P_TO_hybrid = W_hybrid / WP_MN_bat2
P_case1_hy =  W_hybrid/np.array(case1_PW)
P_case2_hy =  W_hybrid/np.array(case2_PW)
P_case3_hy = W_hybrid/np.array(case3_PW)
P_cruise_hy = W_hybrid/np.array(WP_cruise)
S_stall_hybrid = W_hybrid / WS_stall
print('Hybrid P & S Design Parameters: {} ft^2, {} hp'.format(S_hybrid, P_hybrid))

# Hybrid Power vs Wing Area Dimensional Constraint
plt.figure(figsize=(8,4),facecolor='#e3dfd7ff')
ax = plt.gca()  #
ax.set_facecolor("#e3dfd7ff")
plt.axvline(x=S_stall_hybrid, label='Stall', color='k', linestyle='-', linewidth=2)
plt.plot(S_ref_hybrid, np.abs(P_TO_hybrid), label='Takeoff field length', linestyle='-', linewidth=2)
plt.plot(S_ref_hybrid, np.abs(P_MN_hybrid), label="Maneuver", linestyle='-', linewidth=2)
plt.plot(S_ref_hybrid, P_case1_hy,label="Level 1 Climb, CL=1.9, AEO")
plt.plot(S_ref_hybrid, P_case2_hy,label="Level 1 Climb, CL=1.3, OEI")
plt.plot(S_ref_hybrid, P_case3_hy,label="Balked Climb, CL=1.9, AEO")
plt.plot(S_ref_hybrid, P_cruise_hy, label='Cruise', linestyle='-', linewidth=2)
plt.scatter(312, 680, color='red', zorder=5, label='Design for min battery size')
plt.fill_between(S_ref_hybrid, np.maximum(P_TO_hybrid, P_MN_hybrid), 4000, where=(S_ref_hybrid >= S_stall_hybrid),  color='#284e3fff', alpha=0.5, label="Feasible Region")
plt.title('$P - S$ for Hybrid Configuration')
plt.xlabel("$S (ft^2)$")
plt.ylabel("$P (hp)$")
plt.xlim(100,600)
plt.ylim(0,2000)
plt.grid(True)
legend = ax.legend(loc='upper right')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#e3dfd7ff')
plt.show()