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
cargo = 2000#lbs (liquid or solid materials)
passenger = 160#lbs (single person avg)
runway = 1000#ft (runway length)
trapped_fuel_factor = 1.06 #constant multiplied for trapped fuel
composite_factor = 0.95 # empty weight fraction factor if composites are used

# Performance Goals
V_MO = 250#kts = nmi/h (max oper limit speed/max struct cruise speed)
MSS = 100#kts = nmi/h (max stall speed)
MTOW = 19,000#lb (max gross weight/max takeoff weight)

# Conversions
ft2nmi = 1/6076.12 #nmi
acre2nmisqr = 1/847.5 #nmi^2
knots2fps = ((1/ft2nmi)/(60*60)) #ft/s

#Step 2: Obtain Historical/Market Data

# Constants from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
def prop_cbhp(type, #fixed pitch (fp), variable pitch (vp), turboprop (tp)
                condition #cruise or loiter
                ):
    #cruise
    if condition == "cruise" or "Cruise":
        if type == "fp":
            cbhp = 0.4# Piston Prop (fixed pitch)
            np = 0.8# efficiency
            print('Fixed Pitch (Piston)')

        if type == "vp":
            cbhp = 0.4# Piston Prop (variable pitch)
            np = 0.8# efficiency
            print('Variable Pitch (Piston)')

        if type == "tp":
            cbhp = 0.4# Turbo Prop (fixed pitch)
            np = 0.8# efficiency
            print('Turboprop')

        print('Cruise Condition')

    #loiter
    if condition == "loiter" or "Loiter":
        if type == "fp":
            cbhp = 0.5# Piston Prop (fixed pitch)
            np = 0.7# efficiency
            print('Fixed Pitch (Piston)')

        if type == "vp":
            cbhp = 0.5# Piston Prop (variable pitch)
            np = 0.8# efficiency
            print('Variable Pitch (Piston)')

        if type == "tp":
            cbhp = 0.6# Turbo Prop (fixed pitch)
            np = 0.8# efficiency
            print('Turboprop')

        print('Loiter Condition')

    print('C_bhp =',cbhp)
    print('\u03b7_p =',np)

    return cbhp,np #returns the C_bhp and efficiency estimate based on type of prop and flight condition

def prop_SFC(cbhp, #propeller specific fuel consumption
             np, #efficiency
             V=(V_MO + MSS)/2 #ft/s velocity of aircraft
             ):
    C = cbhp*V/(550*np)
    print('C =',C)
    return C #returns C specific fuel consumption for props nased on C_bhp and prop efficiency

def ag_LD_estimates(b, #ft span
          S_wet, #wetted area
          fixed_lg=True #Whether its fixed landing gear or retractable
          ):

    AR_wet = b**2/S_wet #wetted aspect ratio
    print('Wetted AR =',AR_wet)

    # Curve fit equations from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
    if fixed_lg == True:
        LD_max = -0.767*(AR_wet)**2 + 6.43*(AR_wet) + 2.88 #curve fit estimate
        print('Fixed Landing Gear')
    elif fixed_lg != True:
        LD_max = -0.903*(AR_wet)**2 + 7.55*(AR_wet) + 4.20 #curve fit estimate
        print('Retractable Landing Gear')

    # Standard LD equation from Joaquim R. R. A. Martins. The Metabook of Aircraft Design. Feb 2023
    LD = 0.943*LD_max 
    print('L/D =',LD)
    print('L/D_max =',LD_max)
    return LD_max,LD #returns standard L/D and L/D max

def turn(V=(V_MO + MSS)/2, #knts velocity during aerial application is assumed to be average speed (max+min)/2
                g=32.17,#ft/s^2 gravity
                beta=45#deg bank angle
                ):
    V_apply = V*knots2fps #knots --> ft/s conversion
    turning_radius = V_apply^2/(g*np.tan(np.rad2deg(beta))) #ft
    turning_radius_nmi = turning_radius*ft2nmi #ft --> nmi conversion
    turn = np.pi*turning_radius #distance flown for a turn in ft (semi-circle)
    turn_nmi = turn*ft2nmi #ft --> nmi conversion
    print('Bank Angle =',beta,'\u00b0')
    print('Turning Radius =',turning_radius,'ft')
    print('Turning Distance =',turn,'ft')

    return turn_nmi,turning_radius_nmi #returns turning parameters

# Step 3: List Assumptions
# Assumptions
    # application area is a square
    # turning radius taken at average speed = (MSS+V_MO)/2
    # aerial application uses same fuel as cruise
    # assume cargo weight is constant (does not lose weight during flight)

def cruise_range(mission_profile #ferry or standard
          ):
    if mission_profile == "ferry" or "Ferry":
        range =  ferry_range 
        print('Ferry Range =',range,'nmi')
    elif mission_profile == 'cruise' or 'Cruise': # range must be calculated for standard mission profile
        turn_d,turn_rad = turn() #nmi
        des_length = (des_area*acre2nmisqr)**(1/2) #design area length in nmi
        num_of_turns = int(des_length/turn_rad) #number of turns needed to cover area
        range = num_of_turns*(turn_d + des_length) + des_rad*2 #total distance flown in nmi
        print('Standard Range =',range,'nmi')
    else:
        return print('Error, please use Cruise or Ferry')

    return range #returns the total range flown during cruise in nmi

# Step 4: Compute Payload Weight (humans + cargo)

payload = cargo + passenger #lbs

# Step 5: Estimate Fuel Fractions for Cruise & Loiter

def cruise_ff(mission_profile, #ferry or standard
                         propeller_type, #fp, vp, or tp (fixed/variable pitch, or turboprop)
                         b, #span
                         S_wet, #wetted area
                         V=(MSS+V_MO)/2 #knots, cruise speed, default is (Vmin+Vmax)/2
                         ):
    R = cruise_range(mission_profile)
    LD_max, LD = ag_LD_estimates(b,S_wet)
    c_bhp,np = prop_cbhp(propeller_type, condition='cruise')
    c = prop_SFC(c_bhp,np,V)

    cruise_fuel_fraction = np.exp(-R*c/(np*LD))
    print('Cruise Fuel Fraction =',cruise_fuel_fraction)

    return cruise_fuel_fraction
    
def loiter_ff(mission_profile, #ferry or standard
                         propeller_type, #fp, vp, or tp (fixed/variable pitch, or turboprop)
                         b, #span
                         S_wet, #wetted area
                         V=(MSS+V_MO)/2, #knots, cruise speed, default is (Vmin+Vmax)/2
                         E=30#m loiter time
                         ):
    R = cruise_range(mission_profile)
    LD_max, LD = ag_LD_estimates(b,S_wet)
    c_bhp,np = prop_cbhp(propeller_type, condition='cruise')
    c = prop_SFC(c_bhp,np,V)
    E = E/60 #convert to hours

    loiter_fuel_fraction = np.exp(-E*V*c/(np*LD))
    print('Loiter Fuel Fraction =',loiter_fuel_fraction)

    return loiter_fuel_fraction

# Constants from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
def segment_ff(segment, #the segment for the fuel fraction assumption (takeoff, climb, descent, landing)
               ):
    if segment == 'takeoff' or 'Takeoff' or 'Take Off' or 'take off':
        ff = 0.970
    elif segment == 'climb' or 'Climb':
        ff = 0.985
    elif segment == 'descent' or 'Descent':
        ff = 0.990
    elif segment == 'landing' or 'Landing':
        ff = 0.995
    else:
        return print('Segment must be takeoff, climb, descent, or landing')
    
    return ff

# the complete fuel fraction of the mission profile
def complete_ff(mission_profile # ferry or standard
                ):
