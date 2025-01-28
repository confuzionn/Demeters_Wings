# Weight Estimation Code
# Last Updated: 1/22/2025


#-------------------------------------------------------------------------------------------------------------------------#
# Weight Estimation Steps
    # Step 1: Establish Requirements
    # Step 2: Obtain Historical Market Data
    # Step 3: List Assumptions
    # Step 4: Compute Payload Weight (Human + Cargo)


# Step 1: Establish Requirements

# Design Mission Profile
des_rad = 25#nmi (design radius)
des_area = 400#acres (20x20)
reserves = 30#min (additional fuel)
cargo = 2000#lb (liquid or solid materials)
passenger = 160#lb (single person avg)
payload = cargo + passenger

MTOW = 19,000#lb (max gross weight/max takeoff weight)

#Step 2: Obtain Historical/Market Data