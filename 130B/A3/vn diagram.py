# vn diagram code
# B-08-VnDiagram: Slide 12-50


## Design airspeeds ## ---------------------------------------------------------------------------------------------------

# Stalling speed at normal level flight
Vs1 = ((2*W)/rho*S_ref*CL_max)**(1/2)
Vs_1 = ((-2*W)/rho*S_ref*CL_max)**(1/2)

# Design Maneuvering speed or corner speed
    # Speed at which the aircraft will stall at the same point the max_lim load factor is achieved
VA = ?

# Design speed for maximum gust intensity
    # speed at which the aircraft can experience maximum gust loads without structural damage
VB >= Vs1*(1 + ((Kg*Udec*Vc*CLa)/(498*W/S)))

# Positive High AoA (PHA)
    # Lowest speed for n_max (VA)
    # Condition obtained in a pull-up at the highest possible AoA of the wing
    # Resultance force R has a normal component N to the wing and an in-plane C component pointing forward (towards the wind)
    # Maximum C at AoA_max

# Negative High AoA
    # Lowest speed for n_min
    # Condition occurs in intentional maneuvers at low speeds, e.g., entering dive (negative load factor), or due to sudden downdraft at level flight
    # Resultant force R has a normal component N (pointing opposite to PHA) to the wing and in-place C component pointing forward (towards  the wind)
    # Wing assumed to be at negative stalling AOA for steady flow (unlike PHAs)

# Maximum level flgiht speed
    # Speed obtained at maximum continuous thrust in a level, clean (flaps up) configuration
33*(W/S)**(1/2) <= VH <= 0.9*VD

# Design cruising speed
    # Speed selected by the designer to ensure that the aircraft withstands particular loads specific in FARs or applicable AC (i.e., gust loads).
    # *Must be sufficiently greater than VC to account for speed increases due to turbulence
VC >= VB + deltaV

# Design dive speed
    # Maximum speed in a dive at which the structure is design to withstand particular loads specfic in FARs
    # *For transonic aircraft, VD = 1.13*VC
    # *For slower aircraft, VD = 1.25*VC
VD = 1.25*VC

# Positive Low AoA (PLA)
    # Highest speed at n_max (VD)
    # Condition obtained in maximum equivalent airspeed at which the airplane will dieve
    # Resultance force R has a normal component N to the wing and an in-plane C component pointing aft (along with wind)
    # Limit on permissible dive speed depends on type of aircraft, but usually 1.25-1.5 times maximm equivalent speed in level flight

# Negative Low AoA (NLA)
    # Highest speed at n_min
    # Condition occurs in high speed pitch down maneuver or negative gust
    # Resultant force R has a normal component N (pointing opposite to PHA) to the wing and an in-plane C component pointing aft (along the wind)
    # Wing assumed to be at negative AoA limited by highest negative load factor for the aircraft category

## Limit Maneuvering Load Factors ## ---------------------------------------------------------------------------------------------------
    # Positive Limit Maneuvering Load Factor (PLMLF)
        # the PLMLF n for any speed up to VD may not be less than "2.1 + 24,000/(W + 10,000)"...
        # except that n need "not be greater than 3.8" where W is the design maximum taekoff weight

    # Negative Limit Maneuvering Load Factor (NLMLF)
        # The NLMLF may "not be less than 0.4 times the positive load factor" for the normal utility and commute categories
    
    # Maneuvering load factors lower than those specified in this section may be used if...
    # the airplane has design features that make it impossible to exceed these values in flight

## How to choose max limit loads ## ---------------------------------------------------------------------------------------------------
    # This shows the variation in load factor with apirspeed for maneuvers
    # At low speeds the max load factor is constrained by aircraft CL_max ("This is why it is a curve")
n = L/W = (rho_SL*VEAS**2*CL_max)/(2*W/S)

    # At higher speeds it becomes restricted by FAR and AC regulations
    # The maximum maneuver load factor for airplanes weighing more than 50,000 lbs...
    # is usually "n = +2.5". the negative value is "n = -1.0"

# Design Flap Speed
    # May not be less than:
        # 1.6*Vs1 with the flaps in takeoff position at max take off weight (MTOW)
        # 1.8*Vs1 with the flaps in approach position at max landing weight (MLW)
        # 1.8*Vs1 with the flaps in landing position at MLW

## General Comments ## ---------------------------------------------------------------------------------------------------

# Regulations require that strength is shown at those conditons, regardless the feasibility of a particular flight condition
# (i.e., whehter the thrust or power available is equal or greater than the power required;...
# or the flight condition can be achieved by the pilot flying that particular airplane)

# All four conditions must be considered for extreme CG Positions.

# Tail loads must be determined from the most forward and rearward CG positions in the MTOW configuration

# High-lift devices may be critical for wing torsion and shear in the rear spar or downt ail loads -- due to high negative pitching moment

# Note that the maneuvering envelope will not be symmetric because the aircraft are not symmetric (i.e. camber),
# however this is specifically for symmetric loading conditions

# Non-symmetrical loads:
    # non-symmetrical loading conditions (25.347), rolling conditions (25.349), and yawing conditions (25.351) can be important for aerobatics...
    # but less for commercial/commuter planes -- regulations provide fairly conservative design loads;

    # Very important for military fighters -- usually specified by DoD

## Gust Loads ## ---------------------------------------------------------------------------------------------------

# "What if there is a gust?"

# One major concern while considering the loads on your aircraft is the effect of large gusts

# On the edge of storms, aircraft can experience gusts on the range of -1.5g's to 3.5g's

# "How to design safely for that?"

## Angle of Attack is Affected ## ---------------------------------------------------------------------------------------------------

# An instantaneous vertical gust (Ude) changes the AoA of the aircraft:
delta_AoA = arctan(Ude/V) ~= (Ude/V)

# As a result the lift changes:
delta_L = CLa*delta_AoA*q*S = CLa*(Ude/V)*q*S = 1/2*CLa*Ude*rho_0*V*S

# And the load factor:
