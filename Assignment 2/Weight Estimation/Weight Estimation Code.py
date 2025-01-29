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

# Performance Goals
V_MO = 250#kts (max oper limit speed/max struct cruise speed)
MSS = 100#kts (max stall speed)
MTOW = 19,000#lb (max gross weight/max takeoff weight)

#Step 2: Obtain Historical/Market Data
def turning_rad(V=MSS, #ft/s
                g=32.17,#ft/s^2
                beta=45#deg bank angle
                ):
    R = V^2/(g*np.tan(np.rad2deg(beta)))
    return R

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

    return cbhp,np


def prop_SFC(cbhp, #propeller specific fuel consumption
             np, #efficiency
             V #ft/s
             ):
    C = cbhp*V/(550*np)
    return C


