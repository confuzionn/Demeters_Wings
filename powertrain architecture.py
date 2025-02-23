from sympy import symbols, Eq, solve

# POWER BALANCE EQUATIONS
# (1):   P_gt = eff_gt * P_f
# (2):   P_s1 + P_gb = eff_gb * P_gt
# (3):   P_p1 = eff_p1 * P_s1
# (4):   P_e1 = eff_em * P_gb
# (5):   P_e2 = eff_pm * (P_e1 + P_bat)
# (6):   P_s2 = eff_em * P_e2
# (7):   P_p2 = eff_p2 * P_s2
# (8):   shaft_ratio = P_s2 / (P_s1 + P_s2)
# (9):   supplied_ratio = P_bat / (P_bar + P_f)
# (10):  P_p = P_p1 + P_P2

WPp = 10  # Weight to Power Ratio
W = 16856  # Weight

P_p = W / WPp  # Total Propeller Power

eff_GB = 0.96
eff_GT = 0.4
eff_PM = 0.99
eff_EM = 0.96
eff_P = 0.85

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

    solution = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10],
                     (P_gt, P_f, P_s1, P_gb, P_p1, P_e1, P_e2, P_s2, P_p2, P_bat))

    # Print the solution
    print(solution)


get_Powers(0, 1, P_p)
