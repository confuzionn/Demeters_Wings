# Weight Estimation Code
# Last Updated: 1/22/2025

#Standard Imports 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------------------------------------------------#
# Weight Estimation Steps
    # Step 1: Establish Requirements
    # Step 2: Obtain Historical Market Data
    # Step 3: List Assumptions
    # Step 4: Compute Payload Weight (Human + Cargo)


# Step 1: Establish Requirements

# Design Mission Profile
des_rad = 25#nmi (design radius)
ferry_range = 600#nmi
des_area = 400#acres (20x20)
reserves = 30#min (additional fuel)
cargo = 2000#lb (liquid or solid materials)
passenger = 160#lb (single person avg)
runway = 1000#ft (runway length)
payload = cargo + passenger
trapped_fuel_factor = 1.06 #constant multiplied for trapped fuel
composite_factor = 0.95 # empty weight fraction factor if composites are used

# Performance Goals
V_MO = 250#kts = nmi/h (max oper limit speed/max struct cruise speed)
MSS = 100#kts = nmi/h (max stall speed)
MTOW = 19,000#lb (max gross weight/max takeoff weight)

# Conversions
ft2nmi = 1/6076.12 #nmi
acre2nmisqr = 1/847.5 #nmi^2

#Step 2: Obtain Historical/Market Data
def prop_cbhp(type, #fixed pitch (fp), variable pitch (vp), turboprop (tp)
                condition #cruise or loiter
                ):
    #cruise
    if condition == "cruise":
        if type == "fp":
            cbhp = 0.4# Piston Prop (fixed pitch)
            np = 0.8# efficiency

        if type == "vp":
            cbhp = 0.4# Piston Prop (variable pitch)
            np = 0.8# efficiency

        if type == "tp":
            cbhp = 0.4# Turbo Prop (fixed pitch)
            np = 0.8# efficiency

    #loiter
    if condition == "loiter":
        if type == "fp":
            cbhp = 0.5# Piston Prop (fixed pitch)
            np = 0.7# efficiency
        if type == "vp":
            cbhp = 0.5# Piston Prop (variable pitch)
            np = 0.8# efficiency
        if type == "tp":
            cbhp = 0.6# Turbo Prop (fixed pitch)
            np = 0.8# efficiency

    return cbhp,np #returns the C_bhp and efficiency estimate based on type of prop and flight condition


def prop_SFC(cbhp, #propeller specific fuel consumption
             np, #efficiency
             V #ft/s velocity of aircraft
             ):
    C = cbhp*V/(550*np)
    return C #returns C specific fuel consumption for props nased on C_bhp and prop efficiency

def ag_LD_estimates(b, #ft span
          S_wet, #wetted area
          fixed_lg=True #Whether its fixed landing gear or retractable
          ):

    AR_wet = b**2/S_wet #wetted aspect ratio

    # Curve fit equations from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
    if fixed_lg == True:
        LD_max = -0.767*(AR_wet)**2 + 6.43*(AR_wet) + 2.88 #curve fit estimate
    elif fixed_lg != True:
        LD_max = -0.903*(AR_wet)**2 + 7.55*(AR_wet) + 4.20 #curve fit estimate
    
    # Standard LD equation from Joaquim R. R. A. Martins. The Metabook of Aircraft Design. Feb 2023
    LD = 0.943*LD_max 

    return LD_max,LD #returns standard L/D and L/D max

def turn(V_apply=MSS, #ft/s velocity during aerial application
                g=32.17,#ft/s^2 gravity
                beta=45#deg bank angle
                ):
    turning_radius = V_apply^2/(g*np.tan(np.rad2deg(beta)))
    turn = np.pi*turning_radius #distance flown for a turn 
    return turn,turning_radius #returns turning parameters

# Step 3: List Assumptions
# Assumptions
    # application area is a square
    # turning radius taken at max stall speed
    # aerial application uses same fuel as cruise
    # cargo weight is lost during application

def cruise_range(mission_profile #ferry or standard
          ):
    if mission_profile == "ferry":
        range =  ferry_range 
    else: # range must be calculated for standard mission profile
        turn_d,turn_rad = turn()
        des_length = (des_area*acre2nmisqr)**(1/2) #design area length in nmi
        num_of_turns = int(des_length/turn_rad) #number of turns needed to cover area
        range = num_of_turns*(turn_d + des_length) #total distance flown in nmi

    return range #returns the total range flown during cruise in nmi
