import numpy as np
import matplotlib.pyplot as plt

def getZeroLiftDragCoef(WTO, surf):
    c = 1.0447  # From Roskam, for agricultural
    d = 0.5326  # From Roskam, for agricultural

    weight = WTO
    surf_area = surf

    Wing_Loading = weight / surf_area
    print(Wing_Loading)

    # Using Equation log_10(Swet) = c + dlog_10(Weight)
    # It becomes Swet = 10^c + W^d
    Swet = (10 ** c) * (weight ** d)
    print(weight)
    print(f"Wetted Area: {Swet:.2f} ft²")

    # Twin Engines have cf = .0045, while Single Engines have cf = .0055.
    # Since this is a series of smaller engines, value likely falls between twin and single engines

    cf = 0.0050
    f = cf * Swet

    # double check with a and b vals
    # for cf = 0.005, a = -2.301, b = 1
    a = -2.301
    b = 1
    f1 = (10 ** a) * (Swet ** b)

    print(f"Equivalent Parasite Area Comparison: f = {f:.4f}, f₁ = {f1:.4f}")

    CD0 = f / surf_area

    print(f"Zero-Lift Drag Coefficient (CD0): {CD0:.5f}")

    return CD0
def oswald_efficiency(AR, sweep):
    #I am using the equations for estimating Oswald's Efficiency from Raymer
    #The value for straight wing provides a more reasonable approximation, so I am using it for our efficiency value

    sweep_rad = np.radians(sweep)

    e_straight = 1.78 * (1 - 0.045 * AR ** 0.68) - 0.64
    e_swept = 4.61 * (1 - 0.045 * AR ** 0.68) * (np.cos(sweep_rad) ** 0.15) - 3.1

    e = e_straight
    print(f"Oswald's Efficiency Factor (e) for AR = {AR}: {e:.4f}")

    return e

def Cd_Cl_Curve(Cd0, e, AR):
    e_clean = e
    e_takeoff = e - .05
    e_landing = e - .1

    Cd0_clean = Cd0
    Cd0_takeoff = Cd0 + 0.015
    Cd0_takeoff_gear = Cd0_takeoff + 0.02
    Cd0_landing = Cd0 + 0.065
    Cd0_landing_gear = Cd0_landing + 0.02

    #We used a Clmax of 1.3 for cruise, 1.6 for takeoff, and 1.9 for landing
    CL_cruise = np.linspace(-1.3, 1.3, 100)
    CL_takeoff = np.linspace(-1.6, 1.6, 100)
    CL_landing = np.linspace(-1.9, 1.9, 100)

    def CD_clean(CL):
        return Cd0_clean + (CL ** 2) / (np.pi * e_clean * AR)

    def CD_takeoff(CL):
        return Cd0_takeoff + (CL ** 2) / (np.pi * e_takeoff * AR)

    def CD_takeoff_gear(CL):
        return Cd0_takeoff_gear + (CL ** 2) / (np.pi * e_takeoff * AR)

    def CD_landing(CL):
        return Cd0_landing + (CL ** 2) / (np.pi * e_landing * AR)

    def CD_landing_gear(CL):
        return Cd0_landing_gear + (CL ** 2) / (np.pi * e_landing * AR)

    # Calculate Cd values for each condition
    CD_vals_clean = CD_clean(CL_cruise)
    CD_vals_takeoff = CD_takeoff(CL_takeoff)
    CD_vals_landing = CD_landing(CL_landing)
    CD_vals_takeoff_gear = CD_takeoff_gear(CL_takeoff)
    CD_vals_landing_gear = CD_landing_gear(CL_landing)

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(CD_vals_clean, CL_cruise, label="Clean Cruise", color='b')
    plt.plot(CD_vals_takeoff, CL_takeoff, label="Takeoff Flaps, Gear Up", color='g')
    plt.plot(CD_vals_takeoff_gear, CL_takeoff, label="Takeoff Flaps, Gear Down", color='purple')
    plt.plot(CD_vals_landing, CL_landing, label="Landing Flaps, Gear Up", color='r')
    plt.plot(CD_vals_landing_gear, CL_landing, label="Landing Flaps, Gear Down", color='orange')

    plt.xlabel("Drag Coefficient (Cd)")
    plt.ylabel("Lift Coefficient (CL)")
    plt.title("Drag Coefficient vs Lift Coefficient for Different Configurations")
    plt.legend()
    plt.show()


WTO = 16856  # lbs
S = 534.375  # ft^2
AR = 10.03268
sweep = 15 # Degrees

Cd0 = getZeroLiftDragCoef(WTO, S)
e = oswald_efficiency(AR, sweep)

Cd_Cl_Curve(Cd0, e, AR)
