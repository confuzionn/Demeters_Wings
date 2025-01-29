import numpy as np

# ------------------------- Initial Rough Cost Estimate ------------------------- #

def get_rough_aircraft_cost(MTOW):
    b_year = 1989
    t_year = 2035
    b_CEF = 5.17053 + 0.104981 * (b_year - 2025)
    t_CEF = 5.17053 + 0.104981 * (t_year - 2025)
    CEF = t_CEF / b_CEF

    C_aircraft = (10**(1.1846 + 1.2625 * np.log10(MTOW))) * CEF
    return C_aircraft

MTOW = 15000

C_aircraft = get_rough_aircraft_cost(MTOW)
print('Rough Estimated Cost: ${}'.format(C_aircraft))

# ------------------------- RAND DAPCA IV Model for RDT&E and Production Costs ------------------------- #

