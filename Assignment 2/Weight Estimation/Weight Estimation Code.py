# Weight Estimation Code
# Last Updated: 1/22/2025

#Standard Imports 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

#-------------------------------------------------------------------------------------------------------------------------#
# Weight Estimation Steps
    # Step 1: Establish Requirements
    # Step 2: Obtain Historical Market Data
    # Step 3: List Assumptions
    # Step 4: Compute Payload Weight (Human + Cargo)

#-------------------------------------------------------------------------------------------------------------------------#
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
grav = 32.174#ft/s^2

# Conversions
ft2nmi = 1/6076.12 #nmi
acre2nmisqr = 1/847.5 #nmi^2
knots2fps = ((1/ft2nmi)/(60*60)) #ft/s
N2lbf = 0.224809 #lbm*g
m2ft = 3.28084 #ft
kg2lb = 2.20462

#-------------------------------------------------------------------------------------------------------------------------#
#Step 2: Obtain Historical/Market Data

# Constants from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
def prop_cbhp(type=['fp','vp','tp'], #fixed pitch (fp), variable pitch (vp), turboprop (tp)
                condition=['cruise','loiter'] #cruise or loiter
                ):
    
    # print('LOADING: Propeller efficiency...')
    #cruise
    if condition == 'cruise':
        if type == 'fp':
            cbhp = 0.4# Piston Prop (fixed pitch)
            np = 0.8# efficiency

        if type == 'vp':
            cbhp = 0.4# Piston Prop (variable pitch)
            np = 0.8# efficiency

        if type == 'tp':
            cbhp = 0.4# Turbo Prop (fixed pitch)
            np = 0.8# efficiency


    #loiter
    if condition == 'loiter':
        if type == 'fp':
            cbhp = 0.5# Piston Prop (fixed pitch)
            np = 0.7# efficiency

        if type == 'vp':
            cbhp = 0.5# Piston Prop (variable pitch)
            np = 0.8# efficiency

        if type == 'tp':
            cbhp = 0.6# Turbo Prop (fixed pitch)
            np = 0.8# efficiency

    # print('C_bhp =',cbhp)
    # print('\u03b7_p =',np)

    return cbhp,np #returns the C_bhp and efficiency estimate based on type of prop and flight condition

def prop_SFC(cbhp, #propeller specific fuel consumption
             np, #efficiency
             V=(V_MO + MSS)/2 #ft/s velocity of aircraft
             ):
    
    # print('LOADING: Propeller SFC...')
    C = (cbhp*(V*knots2fps)/(550*np))*ft2nmi # 1/nmi
    # print('V =',V,'knots')
    # print('C =',C,'1/nmi')
    return C #returns C specific fuel consumption for props in 1/nmi

def ag_LD_estimates(b, #ft span
          S_wet, #wetted area
          fixed_lg=True #Whether its fixed landing gear or retractable
          ):
    
    # print('LOADING: L/D estimate...')
    # print('Span',b,'ft')
    # print('Wetted Area',S_wet,'ft^2')
    AR_wet = (b**2)/S_wet #wetted aspect ratio
    # print('Wetted AR =',AR_wet)

    # Curve fit equations extrapolated using data from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
    if fixed_lg == True:
        LD_max = 8.59*(AR_wet)**0.56 #curve fit estimate
        # print('Fixed Landing Gear')

    elif fixed_lg == False:
        LD_max = 10.9*(AR_wet)**0.526 #curve fit estimate
        # print('Retractable Landing Gear')

    # Standard LD equation from Joaquim R. R. A. Martins. The Metabook of Aircraft Design. Feb 2023
    LD = 0.943*LD_max 
    # print('L/D =',LD)
    return LD_max,LD #returns standard L/D and L/D max

def turn(V=MSS, #knts velocity during aerial application is assumed to be average speed (max+min)/2
                g=grav,#ft/s^2 gravity
                beta=45#deg bank angle
                ):
    
    # print('LOADING: Turning radius/distance...')
    V_apply = V*knots2fps #knots --> ft/s conversion
    turning_radius = (V_apply**2)/(g*np.tan(np.deg2rad(beta))) #ft
    turning_radius_nmi = turning_radius*ft2nmi #ft --> nmi conversion
    turn = np.pi*turning_radius #distance flown for a turn in ft (semi-circle)
    turn_nmi = turn*ft2nmi #ft --> nmi conversion
    # print('Turning Velocity =',V_apply,'knots')
    # print('Bank Angle =',beta,'\u00b0')
    # print('Turning Radius =',turning_radius,'ft')
    # print('Turning Distance =',turn,'ft')

    return turn_nmi,turning_radius_nmi #returns turning parameters

#-------------------------------------------------------------------------------------------------------------------------#
# Step 3: List Assumptions
# Assumptions
    # application area is a square
    # turning radius taken at average speed = (MSS+V_MO)/2
    # aerial application uses same fuel as cruise
    # assume cargo weight is constant (does not lose weight during flight)

def cruise_range(mission_profile=['ferry','standard'] #ferry or standard
          ):

    # print('LOADING: Cruise range...')
    if mission_profile == 'ferry':
        range =  ferry_range 
        # print('Ferry Range =',range,'nmi')

    elif mission_profile == 'standard': # range must be calculated for standard mission profile
        turn_d,turn_rad = turn() #nmi
        des_length = (des_area*acre2nmisqr)**(1/2) #design area length in nmi
        num_of_turns = int(des_length/turn_rad) #number of turns needed to cover area
        range = num_of_turns*(turn_d + des_length) + des_rad*2 #total distance flown in nmi
        # print('Standard Range =',range,'nmi')
    else:
        return print('Error')

    return range #returns the total range flown during cruise in nmi

#-------------------------------------------------------------------------------------------------------------------------#
# Step 4: Compute Payload Weight (humans + cargo)

payload = cargo + passenger #lbs

#-------------------------------------------------------------------------------------------------------------------------#
# Step 5: Estimate Fuel Fractions for Cruise & Loiter

def cruise_ff(mission_profile=['ferry','standard'], #ferry or standard
                         propeller_type=['fp','vp','tp'], #fp, vp, or tp (fixed/variable pitch, or turboprop)
                         b=any, #span
                         S_wet=any, #wetted area
                         V=(MSS+V_MO)/2, #knots, cruise speed, default is (Vmin+Vmax)/2
                         flg=True #fixed landing gear 
                         ):
    
    # print('LOADING: Cruise fuel fraction...')
    R = cruise_range(mission_profile)
    LD_max, LD = ag_LD_estimates(b=b,S_wet=S_wet,fixed_lg=flg)
    c_bhp,np = prop_cbhp(type=propeller_type, condition='cruise')
    c = prop_SFC(c_bhp,np,V)

    # 
    cruise_fuel_fraction = math.exp(-(R*c/(np*LD)))
    # print('Cruise Fuel Fraction =',cruise_fuel_fraction)

    return cruise_fuel_fraction
    
def loiter_ff(mission_profile=['ferry','standard'], #ferry or standard
                         propeller_type=['fp','vp','tp'], #fp, vp, or tp (fixed/variable pitch, or turboprop)
                         b=any, #span
                         S_wet=any, #wetted area
                         V=MSS, #knots or nmi/h, cruise speed, default is (Vmin+Vmax/2)/2
                         E=30,#m loiter time
                         flg=True #fixed landing gear 
                         ):
    
    # print('LOADING: Loiter range...')
    R = cruise_range(mission_profile)
    LD_max, LD = ag_LD_estimates(b=b,S_wet=S_wet,fixed_lg=flg)
    c_bhp,np = prop_cbhp(type=propeller_type, condition='loiter')
    c = prop_SFC(c_bhp,np,V)
    E = E/60 #convert to hours

    loiter_fuel_fraction = math.exp(-(E*V*c/(np*LD)))
    # print('Loiter Fuel Fraction =',loiter_fuel_fraction)

    return loiter_fuel_fraction

#-------------------------------------------------------------------------------------------------------------------------#
# Step 6: Estimate final fuel fraction
# Constants from Daniel P. Raymer. Aircraft Design: A Conceptual Approach. AIAA, 4th edition, 2006.
def segment_ff(segment=['takeoff','climb','descent','landing'], #the segment for the fuel fraction assumption (takeoff, climb, descent, landing)
               ):
    
    if segment == 'takeoff':
        ff = 0.970
    elif segment == 'climb':
        ff = 0.985
    elif segment == 'descent':
        ff = 0.990
    elif segment == 'landing':
        ff = 0.995
    else:
        return print('Error')
    
    return ff

# the complete fuel fraction (fff) of the mission profile
def complete_ff(mission_profile=['ferry','standard'], #ferry or standard
                propeller_type=['fp','vp','tp'], #fp, vp, or tp (fixed/variable pitch, or turboprop)
                b=any, #span
                S_wet=any, #wetted area
                V=(MSS+V_MO)/2, #knots, cruise speed, default is (Vmin+Vmax)/2
                E=30,#m loiter time
                flg=True #fixed landing gear
                ):
    
    # print('LOADING: Mission profile fuel fraction...')
    if mission_profile == 'ferry': #based on mission profile
        fuel_fract = segment_ff('takeoff') * segment_ff('climb') * \
        cruise_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V,flg=flg) * \
        loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg) * \
        segment_ff('descent') * segment_ff('landing')
        
    elif mission_profile == 'standard': #based on mission profile
        fuel_fract = segment_ff('takeoff') * segment_ff('climb') **2 * \
        cruise_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V,flg=flg) * \
        loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg) * \
        segment_ff('descent') **2 * segment_ff('landing')

    fff = (1-fuel_fract)

    # print('Total Fuel Fraction Wf/W0 =',fff)
    return fff

#-------------------------------------------------------------------------------------------------------------------------#
# Step 7: Iterative Weight Estimation

def weight_estimation(mission_profile=['ferry','standard'], #ferry or standard
                propeller_type=['fp','vp','tp'], #fp, vp, or tp (fixed/variable pitch, or turboprop)
                b=any, #span
                S_wet=any, #wetted area
                V=(MSS+V_MO)/2, #knots, cruise speed, default is (Vmin+Vmax)/2
                E=30, #m loiter time
                A = 0.74, # From Raymer Table 3.1 fuel
                C = -0.03, # From Raymer Table 3.1 fuel
                W0 = 2e4, #lbs, initial empty weight
                energy_density = any, #Wh/kg
                battery_efficiency = 0.7, # Battery efficiency from EAE130A_2025_WQ_Tutorial_Weight_Estimation.html
                flg=True,
                composite=False,
                electric=False,
                name=any
                ):
    print('|-----------------------------------------------------------------------------------------------|')
    print(name)
    # print('LOADING: Weight estimation...')
    W0_history = []   # list of all W0 guesses for plot
    res = 1e-6        # relative convergence tolerance
    error = 2*res     # any value greater than the tolerance
    
    #calculate final fuel fraction
    fff = complete_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg)
    print('Final fuel fraction',fff)

    #find the range traveled in mission profile
    R = cruise_range(mission_profile=mission_profile) 
    
    # Battery specific energy for aviation
    # Source: Barrera, Thomas P., et al. "Next-generation aviation li-ion battery technologies—enabling electrified aircraft." ...
    #         The Electrochemical Society Interface 31.3 (2022): 69.

    # Units: W*h/kg * 3600s/h = Nm/kg --> *lbf/N * ft/m * kg/lb = lbf*ft/lb
    battery_specific_energy = energy_density*3600*N2lbf*m2ft/kg2lb
    print('Battery Specific Energy =',battery_specific_energy,'lbf-ft/lb')
    print('Battery Energy Density =',energy_density/kg2lb,'Wh/lb')
    print('Battery Efficiency = ',battery_efficiency,'%')

    # calculate the LD ratio for electric calculations
    LD_max, LD = ag_LD_estimates(b=b,S_wet=S_wet,fixed_lg=flg)
    AR_wet = b**2/S_wet

    print('Wetted Aspect Ratio =',AR_wet)
    print('L/D =',LD)

    # if there are composites present, the final weight will be lighter
    if composite == True:
        composite = composite_factor
    elif composite == False:
        composite = 1
    else:
        print('Error')
    iter = 0
    if electric == False:
        fff = fff*trapped_fuel_factor
        while error > res:
            iter += 1
            W0_history.append(W0) # add latest value to list
            empty_weight_fraction = A*W0**C*composite #empty weight ratio, lbs/lbs

            # Since the fuel fraction is greater than 1-We/W0, that means I have to consider part of the fuel as the payload
            # converting the fuel into weight, we can then add that to the payload weight
            Wf = fff*W0
            
            W0_new = (payload) / (1 - fff - empty_weight_fraction)
            error = abs(W0_new - W0) / abs(W0_new) #find resolution
            W0 = W0_new #lbs

            # error message
            if iter > 1e4:
                print('ERROR, DID NOT CONVERGE')
                W0 = 0
                empty_weight_fraction = 0
                error = res
            
        W0_history = np.array(W0_history)  # convert list to array
        We = empty_weight_fraction*W0

    else:
        while error > res:
            iter += 1
            W_f = fff*W0
            
            if mission_profile == 'standard':
                Wf_cruise = W_f*cruise_ff(mission_profile='standard', propeller_type=propeller_type,b=b, S_wet=S_wet, V=V,flg=flg)

                # Calculate the weights for the different mission segments
                W1 = segment_ff('takeoff')*W0
                W2 = segment_ff('climb')*W1
                W3 = segment_ff('descent')*W2
                W4 = segment_ff('climb')*W3
                W5 = segment_ff('descent')*W4
                W6 = W5*loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg)

                # Use weights from segments to augment the battery capacity
                E01 = (1-segment_ff('takeoff'))*(W0/Wf_cruise)
                E12 = (1-segment_ff('climb'))*(W1/Wf_cruise)
                E23 = (1-segment_ff('descent'))*(W2/Wf_cruise)
                E34 = (1-segment_ff('climb'))*(W3/Wf_cruise)
                E45 = (1-segment_ff('descent'))*(W4/Wf_cruise)
                E56 = (1-loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg))*(W5/Wf_cruise)
                E67 = (1-segment_ff('landing'))*(W6/Wf_cruise)

                # Use the battery capcity to find a constant to substitute the fuel fraction
                E_constant = (1+E01+E12+E23+E34+E45+E56+E67)

            else:
                Wf_cruise = W_f*cruise_ff(mission_profile='ferry', propeller_type=propeller_type,b=b, S_wet=S_wet, V=V,flg=flg)
                # Calculate the weights for the different mission segments
                W1 = segment_ff('takeoff')*W0
                W2 = segment_ff('climb')*W1
                W3 = segment_ff('descent')*W2
                W4 = W3*loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg)

                # Use weights from segments to augment the battery capacity
                E01 = (1-segment_ff('takeoff'))*(W0/Wf_cruise)
                E12 = (1-segment_ff('climb'))*(W1/Wf_cruise)
                E23 = (1-segment_ff('descent'))*(W2/Wf_cruise)
                E34 = (1-loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg))*(W3/Wf_cruise)
                E45 = (1-segment_ff('landing'))*(W4/Wf_cruise)

                # Use the battery capcity to find a constant to substitute the fuel fraction
                E_constant = (1+E01+E12+E23+E34+E45)

            #Units: lb = ft * lbf * lb/lbf-ft
            m_battery = (R/ft2nmi * W0) / (battery_efficiency*battery_specific_energy*LD) #lbs
            W0_history.append(W0) # add latest value to list                   
            empty_weight_fraction = (A*W0**C)*composite #empty weight ratio, lbs/lbs
            # print('E_constant',E_constant)
            W0_new = (payload) / (1 - empty_weight_fraction - ((m_battery) / W0)*E_constant)   # lbs
            error = abs(W0_new - W0) / abs(W0_new)
            # print(W0)
            W0 = W0_new #lbs 

            # error message
            if iter > 1e4:
                print('ERROR: DID NOT CONVERGE')
                W0 = 0
                empty_weight_fraction = 0
                error = res
            elif W0 < 1:
                print('ERROR: BATTERY ENERGY DENSITY INSUFFICIENT')
                W0 = 0
                empty_weight_fraction = 0
                error = res

        W0_history = np.array(W0_history)  # convert list to array
        We = empty_weight_fraction*W0

        
        W_elec = W0 - payload - We

        print('Battery Weight =',m_battery,'lbs')
        # print('R =',R,'nmi')


    print('Gross Weight =',W0,'lbs')
    print('Empty Weight =',We,'lbs')
    print('|-----------------------------------------------------------------------------------------------|')


    # Plot Convergence
    plt.figure(figsize=(8,4))
    plt.title('{} Weight Estimate Convergence (AR_wet={}, W0={}lbs, We={}lbs)'.format(name,round(AR_wet,1),round(W0),round(We)))
    plt.xlabel('Iteration')
    plt.ylabel('W0 (lbs)')
    plt.plot(W0_history, label='W0', linestyle='-', linewidth=2, marker=None, markersize=8)
    plt.grid(True)
    plt.legend(loc='best')
    plt.show()


def energy_estimation(mission_profile=['ferry','standard'], #ferry or standard
                propeller_type=['fp','vp','tp'], #fp, vp, or tp (fixed/variable pitch, or turboprop)
                b=any, #span
                S_wet=any, #wetted area
                V=(MSS+V_MO)/2, #knots, cruise speed, default is (Vmin+Vmax)/2
                E=30, #m loiter time
                A = 0.74, # From Raymer Table 3.1
                C = -0.03, # From Raymer Table 3.1
                W0 = 2e4, #lbs, initial empty weight
                energy_density = 600, #Wh/kg
                battery_efficiency = 0.85, # Battery efficiency from EAE130A_2025_WQ_Tutorial_Weight_Estimation.html
                flg=True,
                composite=False,
                electric=False,
                name=any
                ):
    print('|-----------------------------------------------------------------------------------------------|')
    print(name)
    # print('LOADING: Weight estimation...')
    W0_history = []   # list of all W0 guesses for plot
    res = 1e-6        # relative convergence tolerance
    error = 2*res     # any value greater than the tolerance
    
    #calculate final fuel fraction
    fff = complete_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg)
    print('Final fuel fraction',fff)

    #find the range traveled in mission profile
    R = cruise_range(mission_profile=mission_profile) 
    
    # Battery specific energy for aviation
    # Source: Barrera, Thomas P., et al. "Next-generation aviation li-ion battery technologies—enabling electrified aircraft." ...
    #         The Electrochemical Society Interface 31.3 (2022): 69.

    # Units: W*h/kg * 3600s/h = Nm/kg --> *lbf/N * ft/m * kg/lb = lbf*ft/lb
    battery_specific_energy = energy_density*3600*N2lbf*m2ft/kg2lb
    print('Battery Specific Energy =',battery_specific_energy,'lbf-ft/lb')
    print('Battery Energy Density =',energy_density/kg2lb,'Wh/lb')
    print('Battery Efficiency = ',battery_efficiency,'%')

    # calculate the LD ratio for electric calculations
    LD_max, LD = ag_LD_estimates(b=b,S_wet=S_wet,fixed_lg=flg)
    AR_wet = b**2/S_wet

    print('Wetted Aspect Ratio =',AR_wet)
    print('L/D =',LD)

    # if there are composites present, the final weight will be lighter
    if composite == True:
        composite = composite_factor
    elif composite == False:
        composite = 1
    else:
        print('Error')
    iter = 0
    if electric == False:
        fff = fff*trapped_fuel_factor
        while error > res:
            iter += 1
            W0_history.append(W0) # add latest value to list
            empty_weight_fraction = A*W0**C*composite #empty weight ratio, lbs/lbs

            # Since the fuel fraction is greater than 1-We/W0, that means I have to consider part of the fuel as the payload
            # converting the fuel into weight, we can then add that to the payload weight
            Wf = fff*W0
            
            W0_new = (payload) / (1 - fff - empty_weight_fraction)
            error = abs(W0_new - W0) / abs(W0_new) #find resolution
            W0 = W0_new #lbs

            # error message
            if iter > 1e4:
                print('ERROR, DID NOT CONVERGE')
                W0 = 0
                empty_weight_fraction = 0
                error = res
            
        W0_history = np.array(W0_history)  # convert list to array
        We = empty_weight_fraction*W0
        print('Final Weight =',W0,'lbs')
        print('|-----------------------------------------------------------------------------------------------|')

        return Wf, W0, We

    else:
        while error > res:
            iter += 1
            W_f = fff*W0
            
            if mission_profile == 'standard':
                Wf_cruise = W_f*cruise_ff(mission_profile='standard', propeller_type=propeller_type,b=b, S_wet=S_wet, V=V,flg=flg)

                # Calculate the weights for the different mission segments
                W1 = segment_ff('takeoff')*W0
                W2 = segment_ff('climb')*W1
                W3 = segment_ff('descent')*W2
                W4 = segment_ff('climb')*W3
                W5 = segment_ff('descent')*W4
                W6 = W5*loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg)

                # Use weights from segments to augment the battery capacity
                E01 = (1-segment_ff('takeoff'))*(W0/Wf_cruise)
                E12 = (1-segment_ff('climb'))*(W1/Wf_cruise)
                E23 = (1-segment_ff('descent'))*(W2/Wf_cruise)
                E34 = (1-segment_ff('climb'))*(W3/Wf_cruise)
                E45 = (1-segment_ff('descent'))*(W4/Wf_cruise)
                E56 = (1-loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg))*(W5/Wf_cruise)
                E67 = (1-segment_ff('landing'))*(W6/Wf_cruise)

                # Use the battery capcity to find a constant to substitute the fuel fraction
                E_constant = (1+E01+E12+E23+E34+E45+E56+E67)

            else:
                Wf_cruise = W_f*cruise_ff(mission_profile='ferry', propeller_type=propeller_type,b=b, S_wet=S_wet, V=V,flg=flg)
                # Calculate the weights for the different mission segments
                W1 = segment_ff('takeoff')*W0
                W2 = segment_ff('climb')*W1
                W3 = segment_ff('descent')*W2
                W4 = W3*loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg)

                # Use weights from segments to augment the battery capacity
                E01 = (1-segment_ff('takeoff'))*(W0/Wf_cruise)
                E12 = (1-segment_ff('climb'))*(W1/Wf_cruise)
                E23 = (1-segment_ff('descent'))*(W2/Wf_cruise)
                E34 = (1-loiter_ff(mission_profile=mission_profile, propeller_type=propeller_type,b=b, S_wet=S_wet, V=V, E=E,flg=flg))*(W3/Wf_cruise)
                E45 = (1-segment_ff('landing'))*(W4/Wf_cruise)

                # Use the battery capcity to find a constant to substitute the fuel fraction
                E_constant = (1+E01+E12+E23+E34+E45)

            #Units: lb = ft * lbf * lb/lbf-ft
            m_battery = (R/ft2nmi * W0) / (battery_efficiency*battery_specific_energy*LD) #lbs
            W0_history.append(W0) # add latest value to list                   
            empty_weight_fraction = (A*W0**C)*composite #empty weight ratio, lbs/lbs
            # print('E_constant',E_constant)
            W0_new = (payload) / (1 - empty_weight_fraction - ((m_battery) / W0)*E_constant)   # lbs
            error = abs(W0_new - W0) / abs(W0_new)
            # print(W0)
            W0 = W0_new #lbs 

            # error message
            if iter > 1e4:
                print('ERROR: DID NOT CONVERGE')
                W0 = 0
                empty_weight_fraction = 0
                error = res
            elif W0 < 1:
                print('ERROR: BATTERY ENERGY DENSITY INSUFFICIENT')
                W0 = 0
                empty_weight_fraction = 0
                error = res

        W0_history = np.array(W0_history)  # convert list to array
        We = empty_weight_fraction*W0
        
        W_elec = W0 - payload - We

        print('Battery Weight =',m_battery,'lbs')
        print('Final Weight =',W0,'lbs')
        print('|-----------------------------------------------------------------------------------------------|')
        # print('R =',R,'nmi')

        return m_battery, W0, We

#-------------------------------------------------------------------------------------------------------------------------#
# Running the code
# weight_estimation('ferry', #ferry or standard mission profile
#                 'vp', #fp, vp, or tp (fixed/variable pitch, or turboprop)
#                 73.5, #span, ft
#                 534.37, #wetted area, ft^2
#                 flg=True, #Landing Gear: True (fixed) or False (retractable)
#                 composite=True, #Whether composites is used or not
#                 electric=True, # Whether its electric or not
#                 energy_density = 550, #Wh/kg, battery energy density
#                 battery_efficiency = 0.85, #factor between 0 and 1
#                 A = 0.459, # From electric aircraft data
#                 C = -0.000268, # From electric aircraft data
#                 name='Fully Electric'
#                 )

def hybrid_estimation(hybrid_ratio= 0.5 # 1 = 100% electric
                      ):
    
    

    m_bat, WO_e, WE_e = energy_estimation('ferry', #ferry or standard mission profile
                'vp', #fp, vp, or tp (fixed/variable pitch, or turboprop)
                73.5, #span, ft
                534.37, #wetted area, ft^2
                flg=True, #Landing Gear: True (fixed) or False (retractable)
                composite=True, #Whether composites is used or not
                electric=True, # Whether its electric or not
                energy_density = 650, #Wh/kg, battery energy density
                battery_efficiency = 0.8, #factor between 0 and 1
                A = 0.459, # From electric aircraft data
                C = -0.000268, # From electric aircraft data
                name='Electric'
                )
    
    m_fuel, WO_f, WE_f = energy_estimation('ferry', #ferry or standard mission profile
                'vp', #fp, vp, or tp (fixed/variable pitch, or turboprop)
                73.5, #span, ft
                534.37, #wetted area, ft^2
                flg=True, #Landing Gear: True (fixed) or False (retractable)
                composite=True, #Whether composites is used or not
                electric=False,
                name='Fuel'
                )
    
    # averaging the empty weight between two energy modes
    WE = hybrid_ratio*WE_e + (1-hybrid_ratio)*WE_f
    Energy_weight = hybrid_ratio*m_bat + (1-hybrid_ratio)*m_fuel
    W0 = WE + Energy_weight + payload

    print('|-----------------------------------------------------------------------------------------------|')
    print('Hybrid')
    print('Hybrid Fraction = {}% Electric'.format(hybrid_ratio*100))
    print('Empty Weight =',WE,'lbs')
    print('Battery Weight =',hybrid_ratio*m_bat,'lbs')
    print('Fuel Weight =',(1-hybrid_ratio)*m_fuel,'lbs')
    print('Total Weight =',W0,'lbs')
    print('|-----------------------------------------------------------------------------------------------|')

hybrid_estimation(hybrid_ratio= 0.5 # 1 = 100% electric
                      )
