from Code.Codebase import *

def Use_Raw(**kwargs): #Load the Raw Run Commands
    ad = ADaX(planename=None) #Create object just for using run()
    raw_file_contents = open("Raw Run Commands",'r').read() #Open Raw Run Commands and get string of contents
    ad.input(raw_file_contents) #Put the string contents into the input list
    ad.run(print_output=True, Log=True, autodebug=True) #Run it
    return None

def Trim(ad: ADaX, **kwargs):
    try: ad.flags['output_configured'] #See if output has been configured
    except: ad.output_config(['t Alpha', 't e', 't CLtot', 't CDind', 's Xnp', 't Elevator']) #If not, default output
    ad._load() # load
    ad._oper() # enter oper
    ad.vcv('d1 pm 0') # trim vcv
    ad.vcv('a c 0.43') # additional stuff
    ad._x() # execute flow calculation
    ad._save() # save outputs
    ad.run(print_output=True, Log=True, autodebug=True) # run it
    return None

def Show(ad: ADaX, **kwargs): #⚠⚠⚠ will timeout AVL, programmed to exit after 60 seconds. use CTRL+C in Visual Studio to close manually.
    ad._load() # load plane
    ad._geometry() # enter oper, open geometry, enter keystroke menu (AVL dies if you do this idk why)
    ad.run(print_output=False, Log=False, autodebug=False) # runs without printing output, writing log, or autodebugging
    return None

def Stall_Prediction(ad: ADaX, Clx: float, **kwargs):
    """
    \nStall prediction based on airfoil maximum Cl (Clx)
    \nChecks strip forces along wing
    \nMakes the assumption that if strip Cl reaches 95% of Clx, wing stalls
    \n⚠ AVL is a potential flow solver; it cannot accurately predict stall! ⚠
    """
    Clx = 0.95*Clx # use 95% of airfoil Clx
    a = 0 # angle of attack
    k = -10 # fixed point iteration constant
    Cl_error = -1 # placeholder error value
    print("Beginning stall prediction iteration with Clx = {:.3} (95% of {:.3})".format(Clx, Clx/0.95)) #printing
    header = ['iter no.', 'Alpha', 'Max strip Cl', 'res']
    print('{:>10}  {:>10}  {:>13}  {:>10}'.format(*header))
    print('─'*55)
    n = 0
    while abs(Cl_error) > 1e-3: # fixed point iteration
        a += Cl_error*k # update alpha
        n += 1 # increment iteration number
        ad._load() # running avl at alpha
        ad._oper()
        ad.vcv('a a {}'.format(a))
        ad._x()
        ad._save()
        ad.run(print_output=False,Log=False,autodebug=False)
        content = open('AVL/Output/Strip forces').readlines() # read the strip forces file
        cl = []
        for i, line in enumerate(content): # basically, looking for the highest Cl on the wing.
            if ' j ' in line:
                for j,line in enumerate(content[i+1:]):
                    if line == '\n': 
                        break
                    else: 
                        cl.append(line.split()[9])
                cl = max(np.array(list(map(float,cl))))
                break
        Cl_error = cl - Clx # error is difference between airfoil max and maximum found on wing
        header = [n, a, cl, Cl_error]
        print('{:>10}  {:>10.4f}  {:>13.4f}  {:>10.5f}'.format(*header))
    print("Finished stall prediction at Alpha = {}".format(np.around(a,4)))
    return None

def Inviscid_Polar(ad: ADaX, a_array: np.ndarray, trim_vcv: str, **kwargs):
    """
    \nGenerates a drag polar file for a range of alpha
    \n⚠ Name your elevator "Elevator" for best results
    \nTakes a_array as a numpy array of desired angles of attack
    \nTakes a trim vcv, should be 'd3 pm 0' or something, where d3 is the elevator
    \nDoes not return anything, generates "AVL/Output/Invis Polar.dat"
    \n
    \nRoom for improvement: write results to "Invis Polar.dat" after each run;
    \ncurrently it attempts to write it all at the end. If there's any error,
    \nnothing gets written.
    """
    ad.output_config(['t Alpha', 't CLtot', 't CDind', 't Elevator']) # specific output configuration
    polar_dict = ad.sweep(vc='a a', v_array=a_array, executes=[('vcv', trim_vcv)]) # gets a dict[dict] from the sweep()
    f = open('AVL/Output/Invis Polar.dat', 'w') # open a file to write results
    f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format('Alpha', 'CLtot', 'CDind', 'Elevator'))
    f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format('-'*20, '-'*20, '-'*20, '-'*20))
    for k, val in polar_dict.items():
        f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format(val['Alpha'], val['CLtot'], val['CDind'], val['Elevator']))
    return None