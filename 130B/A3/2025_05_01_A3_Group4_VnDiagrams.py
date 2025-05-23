# vn diagram code
# B-08-VnDiagram: Slide 12-50

# imports
import numpy as np
import matplotlib.pyplot as plt

# ## Known Variables ## ---------------------------------------------------------------------------------------------------
CLmax = 1.9
CLmin = -0.8
# MTOW = 8145 #lbs
MTOW = 5507+160+95+160 #lbs
S_ref = 312 #ft^2
c_ref = 6 #ft
altitude = 8000 #ft
rho_8k = 0.00187 #slugs/ft^3
rho_sl = 0.0023769 #slugs/ft^3
grav = 32.17 #ft/s^2
n_max = 4.4
n_min = -1.8
VEAS_max = 250 #knots
VEAS_min = 0 #knots
e_efficiency = 0.7922
Aspect_ratio = 8.7
iterations = 100

# ## Design airspeeds ## ---------------------------------------------------------------------------------------------------

# # Stalling speed at normal level flight
# Vs1 = ((2*W)/rho*S_ref*CL_max)**(1/2)
# Vs_1 = ((-2*W)/rho*S_ref*CL_max)**(1/2)

# DEFINING STALLING SPEED EQUATION
def VS1(
    n, #load factor
    W, #lbs
    S, #ft^2
    rho, #slug/ft^3
    CL_max
    ):
    Vs1 = np.sqrt((2*n*W)/(rho*S*CL_max))

    return Vs1

def VS_1(
    n, #load factor
    W, #lbs
    S, #ft^2
    rho, #slug/ft^3
    CL_max
    ):
    Vs_1 = np.sqrt((2*-n*W)/(rho*S*CL_max))

    return Vs_1

# # Design Maneuvering speed or corner speed
#     # Speed at which the aircraft will stall at the same point the max_lim load factor is achieved
# VA = ?
def VA(
    n,
    CL_max, 
    rho_SL,  #slug/ft^3
    W, #lbs
    S   #ft^2
    ):
    n_max = 4.4
    V_A = np.sqrt(n_max * (2 * (n*W)/ S) / (rho_SL * CL_max))

    return V_A

def VG(
    n,
    CL_min, 
    rho_SL,  #slug/ft^3
    W, #lbs
    S   #ft^2
    ):
    n_min = -1.8
    V_G = np.sqrt(-n_min * (2 * (n*W) / S) / (rho_SL * CL_min))

    return V_G

# # Design speed for maximum gust intensity
#     # speed at which the aircraft can experience maximum gust loads without structural damage
# VB >= Vs1*(1 + ((Kg*Ude*Vc*CLa)/(498*W/S)))

# DEFINING DESIGN SPEED FOR MAX GUST EQUATION
def VB(
    n, 
    Vs1, #knots
    Kg,
    Udec, #ft/s
    Vc, #knots
    CLa, #rad^-1
    W, #lbs
    S, #ft^2
    ):
    # Greater than or equal to
    V_B = Vs1*(1 + ((Kg*Udec*Vc*CLa)/(498*(n*W)/S)))
    return V_B

# # Positive High AoA (PHA)
#     # Lowest speed for n_max (VA)
#     # Condition obtained in a pull-up at the highest possible AoA of the wing
#     # Resultance force R has a normal component N to the wing and an in-plane C component pointing forward (towards the wind)
#     # Maximum C at AoA_max

# # Negative High AoA
#     # Lowest speed for n_min
#     # Condition occurs in intentional maneuvers at low speeds, e.g., entering dive (negative load factor), or due to sudden downdraft at level flight
#     # Resultant force R has a normal component N (pointing opposite to PHA) to the wing and in-place C component pointing forward (towards  the wind)
#     # Wing assumed to be at negative stalling AOA for steady flow (unlike PHAs)

# # Maximum level flight speed
#     # Speed obtained at maximum continuous thrust in a level, clean (flaps up) configuration
# 33*(W/S)**(1/2) <= VH <= 0.9*VD

# DEFINING MAX LEVEL FLIGHT SPEED EQUATION
def VH(
    n, 
    W, #lbs
    S, #ft^2
    ):
    # Greater than or equal to
    V_H = 33*np.sqrt((n*W)/S)
    return V_H

def VD(
    V_H #knots
    ):
    # Greater than or equal to
    V_D = V_H/(0.9)
    return V_D



# # Design cruising speed
#     # Speed selected by the designer to ensure that the aircraft withstands particular loads specific in FARs or applicable AC (i.e., gust loads).
#     # *Must be sufficiently greater than VC to account for speed increases due to turbulence
# VC >= VB + deltaV

# # Design dive speed
#     # Maximum speed in a dive at which the structure is design to withstand particular loads specfic in FARs
#     # *For transonic aircraft, VD = 1.13*VC
#     # *For slower aircraft, VD = 1.25*VC
# VD = 1.25*VC

# DEFINING CRUISING SPEED EQUATION
def VC(
    VD #knots
    ):
    V_C = VD/1.25
    return V_C

# # Positive Low AoA (PLA)
#     # Highest speed at n_max (VD)
#     # Condition obtained in maximum equivalent airspeed at which the airplane will dieve
#     # Resultance force R has a normal component N to the wing and an in-plane C component pointing aft (along with wind)
#     # Limit on permissible dive speed depends on type of aircraft, but usually 1.25-1.5 times maximum equivalent speed in level flight

# # Negative Low AoA (NLA)
#     # Highest speed at n_min
#     # Condition occurs in high speed pitch down maneuver or negative gust
#     # Resultant force R has a normal component N (pointing opposite to PHA) to the wing and an in-plane C component pointing aft (along the wind)
#     # Wing assumed to be at negative AoA limited by highest negative load factor for the aircraft category

# ## Limit Maneuvering Load Factors ## ---------------------------------------------------------------------------------------------------
#     # Positive Limit Maneuvering Load Factor (PLMLF)
#         # the PLMLF n for any speed up to VD may not be less than "2.1 + 24,000/(W + 10,000)"...
#         # except that n need "not be greater than 3.8" where W is the design maximum taekoff weight

#     # Negative Limit Maneuvering Load Factor (NLMLF)
#         # The NLMLF may "not be less than 0.4 times the positive load factor" for the normal utility and commute categories
    
#     # Maneuvering load factors lower than those specified in this section may be used if...
#     # the airplane has design features that make it impossible to exceed these values in flight

# # Design Flap Speed
#     # VF may not be less than:
#    (1)  # 1.6*Vs1 with the flaps in takeoff position at max take off weight (MTOW)
#    (2)  # 1.8*Vs1 with the flaps in approach position at max landing weight (MLW)
#    (3)  # 1.8*Vs1 with the flaps in landing position at MLW

# DEFINING DESIGN FLAP SPEED EQUATION
def VF(
    Vs1, #knots
    case_num #look comments above
    ):
    if case_num == 1:
        V_F = 1.6*Vs1
    elif case_num == 2:
        V_F = 1.8*Vs1
    elif case_num == 3:
        V_F = 1.8*Vs1
    else:
        print('Not Valid Case Number')
    
    return V_F

# ## General Comments ## ---------------------------------------------------------------------------------------------------

# # Regulations require that strength is shown at those conditons, regardless the feasibility of a particular flight condition
# # (i.e., whehter the thrust or power available is equal or greater than the power required;...
# # or the flight condition can be achieved by the pilot flying that particular airplane)

# # All four conditions must be considered for extreme CG Positions.

# # Tail loads must be determined from the most forward and rearward CG positions in the MTOW configuration

# # High-lift devices may be critical for wing torsion and shear in the rear spar or downt ail loads -- due to high negative pitching moment

# # Note that the maneuvering envelope will not be symmetric because the aircraft are not symmetric (i.e. camber),
# # however this is specifically for symmetric loading conditions

# # Non-symmetrical loads:
#     # non-symmetrical loading conditions (25.347), rolling conditions (25.349), and yawing conditions (25.351) can be important for aerobatics...
#     # but less for commercial/commuter planes -- regulations provide fairly conservative design loads;

#     # Very important for military fighters -- usually specified by DoD

# ## Gust Loads ## ---------------------------------------------------------------------------------------------------

# # "What if there is a gust?"

# # One major concern while considering the loads on your aircraft is the effect of large gusts

# # On the edge of storms, aircraft can experience gusts on the range of -1.5g's to 3.5g's

# # "How to design safely for that?"

# ## How to choose max limit loads ## ---------------------------------------------------------------------------------------------------
#     # This shows the variation in load factor with airspeed for maneuvers
#     # At low speeds the max load factor is constrained by aircraft CL_max ("This is why it is a curve")
# n = L/W = (rho_SL*VEAS**2*CL_max)/(2*W/S)

# DEFINING MAXIMUM LIMIT LOAD EQUATION
def max_lim_load(
    rho_SL, #slugs/ft^3
    VEAS, #knots
    CL_max,
    W, #lbs
    S, #ft^2
    ):
    n_maxlim = (rho_SL*VEAS**2*CL_max)/(2*W/S)
    return n_maxlim

def min_lim_load(
    rho_SL, #slugs/ft^3
    VEAS, #knots
    CL_max,
    W, #lbs
    S, #ft^2
    ):
    n_minlim = (rho_SL*VEAS**2*CL_max)/(2*-W/S)
    return n_minlim

#     # At higher speeds it becomes restricted by FAR and AC regulations
#     # The maximum maneuver load factor for airplanes weighing more than 50,000 lbs...
#     # is usually "n = +2.5". the negative value is "n = -1.0"

# ## Angle of Attack is Affected ## ---------------------------------------------------------------------------------------------------

# # An instantaneous vertical gust (Ude) changes the AoA of the aircraft:
# delta_AoA = arctan(Ude/V) ~= (Ude/V) 

# # As a result the lift changes:
# delta_L = CLa*delta_AoA*q*S = CLa*(Ude/V)*q*S = 1/2*CLa*Ude*rho_0*V*S

# # And the load factor:
# n = 1 + delta_n = 1 + delta_L/W = 1 + (CLa*Ude*rho_0*V)/(2*W/S)

# DEFINING INSTANTANEOUS LOAD FACTOR CHANGE DUE TO GUST EQUATION
def n_insta_gust(
    CLa, #rad^-1
    Ude, #ft/s
    rho_0, #slug/ft^3
    V, #knots
    W, #lbs
    S #ft^2
    ):
    n_instagust = 1 + (CLa*Ude*rho_0*V)/(2*W/S)
    return n_instagust

# # To account for non-instantaneous nature of gusts we add a factor Kg
# # The gust load may be computed from the expresssion given in FAR part 25 (Part 23 uses the same)
# n = 1 +- (Kg*CLa*Ude*VEAS)/(498*W/S)

# DEFINING NON-INSTANTANEOUS LOAD FACTOR CHANGE DUE TO GUST EQUATION
def n_noinsta_gust_pos(
    Kg,
    CLa, #rad^-1
    Ude, #ft/s
    VEAS, #knots
    W, #lbs
    S #ft^2
    ):
    n_gustpos = 1 + (Kg*CLa*Ude*VEAS)/(498*W/S)
    return n_gustpos

def n_noinsta_gust_neg(
    Kg,
    CLa, #rad^-1
    Ude, #ft/s
    VEAS, #knots
    W, #lbs
    S #ft^2
    ):
    n_gustneg = 1 - (Kg*CLa*Ude*VEAS)/(498*W/S)
    return n_gustneg

# # Definitions
# # Ude = Reference gust, EAS (ft/s)
# # VEAS = Equivalent airspeed (knots)
# # CLa = Lift slope in rad^-1
# # Kg = Gust alleviation factor
# Kg = (0.88*mu)/(5.3 + mu)
# mu = (2*W/S)/(rho*c_bar*CLa*g)

# DEFINING REFERENCE GUST EQUATION
def UDE_VC(
    h, #ft
    ):
    Ude_for_VC = 56 + ((44 - 56)/15000)*h
    return Ude_for_VC

# DEFINING KG EQUATION
def KG(
    W, #lbs
    S, #ft^2
    rho, #slug/ft^3
    c_bar, #ft
    CLa, #rad^-1
    g #ft/s^2
    ):
    mu = (2*W/S)/(rho*c_bar*CLa*g)
    Kg = (0.88*mu)/(5.3 + mu)
    return Kg

#  DEFINING CLA EQUATION

def CLa(
    AR, 
    e
    ):
    C_la = (2 * np.pi) * (np.pi / 180)
    C_La = C_la / (1 + 57.3 * C_la / (np.pi * AR * e))
    return C_La

# CALCULATING CONSTANTS
CL_a = CLa(Aspect_ratio,e_efficiency)
Kg = KG(MTOW,S_ref,rho_8k,c_ref,CL_a,grav)
Ude_VC = UDE_VC(altitude)
Ude_VB = Ude_VC
Ude_VD = Ude_VC/2

# CALCULATING VELOCITIES BASED ON LOAD FACTOR
n = np.linspace(n_min,n_max,iterations+1)

# Velocity history
Vs1 = []
Vs_1 = []
V_H = []
V_D = []
V_C = []
V_B = []
V_F_case1 = []
V_F_case2 = []
V_F_case3 = []
V_A = []
V_G = []

count = -1
for i in n:
    count = count + 1
    Vs1.append(VS1(i,MTOW,S_ref,rho_8k,CLmax))
    Vs_1.append(VS_1(i,MTOW,S_ref,rho_8k,CLmax))
    V_H.append(VH(i,MTOW,S_ref))
    V_D.append(VD(V_H[count]))
    V_C.append(VC(V_D[count]))
    V_B.append(VB(i,Vs1[count],Kg,Ude_VB,V_C[count],CL_a,MTOW,S_ref))
    V_F_case1.append(VF(Vs1[count],1))
    V_F_case2.append(VF(Vs1[count],2))
    V_F_case3.append(VF(Vs1[count],3))
    V_A.append(VA(i,CLmax,rho_sl,MTOW,S_ref))
    V_G.append(VG(i,CLmin,rho_sl,MTOW,S_ref))

# CALCULATING LOAD FACTOR BASED ON VELOCITIES
V_EAS = np.linspace(VEAS_min,VEAS_max,iterations+1)

# Load Factor history
n_maxlim = []
n_minlim = []
n_instagust = []
n_gustpos = []
n_gustneg = []

count = -1
for v in V_EAS:
    count = count + 1
    n_maxlim.append(max_lim_load(rho_sl,v,CLmax,MTOW,S_ref))
    n_minlim.append(min_lim_load(rho_sl,v,CLmax,MTOW,S_ref))
    n_instagust.append(n_insta_gust(CL_a,Ude_VC,rho_sl,v,MTOW,S_ref))
    n_gustpos.append(n_noinsta_gust_pos(Kg,CL_a,Ude_VC,v,MTOW,S_ref))
    n_gustneg.append(n_noinsta_gust_neg(Kg,CL_a,Ude_VC,v,MTOW,S_ref))

# PLOT THE V-N DIAGRAM
plt.figure()
plt.xlabel('Velocity [knots]')
plt.ylabel('Load Factor, n')
plt.grid()
ax = plt.gca()
ax.set_xlim([VEAS_min, VEAS_max])
ax.set_ylim([n_min, n_max])
plt.title(['V-N Diagram: Weight =',MTOW,'lbs'])

# Plotting calculated V
plt.plot(Vs1,n,'b:',label='VS1')
plt.plot(Vs_1,n,'b:',label='VS-1')
plt.plot(V_H,n,':',label='VH')
plt.plot(V_D,n,':',label='VD')
plt.plot(V_C,n,':',label='VC')
plt.plot(V_B,n,':',label='VB')
plt.plot(V_F_case1,n,'g-.',label='VF1')
plt.plot(V_F_case2,n,'g-',label='VF2')
plt.plot(V_F_case3,n,'g:',label='VF3')
plt.plot(V_A,n,':',label='VA')
plt.plot(V_G,n,':',label='VG')

# plotting calculated n
plt.plot(V_EAS,n_maxlim,'k--',label='Max LF')
plt.plot(V_EAS,n_minlim,'k--',label='Min LF')
plt.plot(V_EAS,n_instagust,'m--',label='Instant Gust LF')
plt.plot(V_EAS,n_gustpos,'r--',label='Gust +LF')
plt.plot(V_EAS,n_gustneg,'r--',label='Gust -LF')

plt.legend()

plt.show()