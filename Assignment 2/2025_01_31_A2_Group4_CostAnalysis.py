import numpy as np

# ------------------------- Initial Rough Cost Estimate ------------------------- #

def get_rough_aircraft_cost(MTOW, SHP_TO):
    b_year = 1989
    t_year = 2024
    b_CEF = 5.17053 + 0.104981 * (b_year - 2024)
    t_CEF = 5.17053 + 0.104981 * (t_year - 2024)
    CEF = t_CEF / b_CEF

    C_aircraft = (10**(1.1846 + 1.2625 * np.log10(MTOW))) * (CEF)
    C_engines = 10**(2.5262 + 0.9465 * np.log10(SHP_TO)) * (CEF)
    C_airframe = C_aircraft - C_engines
    return C_aircraft, C_engines, C_airframe, CEF

C_aircraft, C_engine, C_airframe, CEF = get_rough_aircraft_cost(15000, 1300)
print('Rough Estimated Cost: ${}'.format(C_aircraft))

# ------------------------- RAND DAPCA IV Model for RDT&E and Production (10% Flyaway) Costs ------------------------- #

def RAND_DAPCAP_IV(W_e, V, Q):
    HE = 4.86 * (W_e**0.777) * (V**0.894) * (Q**0.163)      # Engineering hours
    HT = 5.99 * (W_e**0.777) * (V**0.696) * (Q**0.263)      # Tooling hours
    HM = 7.37 * (W_e**0.82) * (V**0.0484) * (Q**0.641)      # Manufacturing hours
    HQ = 0.133 * HM                                         # Quality control hours
    CD = 91.3 * (W_e**0.630) * (V**1.3)                     # Development support cost
    FTA = 2                                                 # Number of flight test aircraft
    CF = 2498 * (W_e**0.325) * (V**0.822) * (FTA**1.21)     # Flight test cost
    CM = 22.1 * (W_e**0.921) * (V**0.621) * (Q**0.799)      # Manufacturing materials cost
    
    # Hourly rate linear fit equations
    y = 2024              # Production year
    RT = 2.883 * y - 5666 # Tooling rate
    RE = 2.576 * y - 5058 # Engineering rate
    RQ = 2.60 * y - 5112  # Quaility rate
    RM = 2.316 * y - 4552 # Manufacturing rate

    RDTE = (HE * RE) + (HT * RT) + (HM * RM) + (HQ * RQ) + CD + CF + CM
    #flyaway = (C_eng * N_eng) + C_avionics
    return RDTE
RDTE_flyaway = RAND_DAPCAP_IV(10000, 250, 10)
print(RDTE_flyaway)

# ------------------------- Operation and Maintenance Costs (DOC) ------------------------- #

def operation_maintenance(t_b, CEF, MTOW, SHP_TO, W_A, R, RL):
    C_crew = (440 + 0.532 * (MTOW/1000) * (CEF) * (t_b))
    C_airport = 1.5 * (MTOW/1000) * (CEF)
    C_navigation = 0.5 * (CEF) * (1.852*R)/t_b * np.sqrt((0.00045359237 * MTOW)/50)
    C_ML = 1.03 * (3 + (0.067 * W_A)/1000) * RL
    C_MM = 1.03 * (30 * CEF) + (0.79 * 10**-5 * C_airframe)
    C_airframe_maintenance = (C_ML + C_MM) * t_b
    C_EM = 1.03 * 1.3 * (0.4956 + 0.0532 * (SHP_TO)/1000 * 1100/5000 + 0.1) * RL
    
    U_annual = 1.5 * 10**3 * (3.4546 * t_b + 2.994 - (12.289 * t_b**2 - 5.6626 * t_b + 8.964)**0.5)
    C_insurance = ((0.02 * C_aircraft)/U_annual) * t_b
    DOC = C_crew + C_airport + C_navigation + C_ML + C_MM + C_airframe_maintenance + C_EM + C_insurance
    C_registration = (0.001 + 10**-8 * MTOW) * DOC
    C_finance = 0.07 * DOC
    totalDOC = DOC + C_registration + C_finance

    return totalDOC

totalDOC = operation_maintenance(4, CEF, 17000, 1300, 7000, 600, 100)
print(totalDOC)