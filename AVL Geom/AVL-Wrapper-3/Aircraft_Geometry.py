from Code.Codebase import *

def Update_Geometry(**kwargs):
    """
    Rewrites your .avl file when called.
    Any of this code may be changed to customize how the .avl file is generated, based on whichever parameters you choose.
    If you want to know how this stuff is organized, look into Documentation/Geometry.ipynb
    """
    #Planform parameters for defining your main wing
    S = 8890
    ar = 6.89
    tr = 0.4
    #Four angles for defining your wing
    le_sweep = 30 /180*np.pi
    dihedral = -3 /180*np.pi
    incidence = 0 
    twist = -3
    #Calculating values more relevant to AVL
    b = np.sqrt(S*ar)
    cr = 2*S/(b*(1+tr))
    ct = tr*cr
    #Airfoil files
    af1 = 'AVL/Airfoils/sc20710.dat'
    af2 = 'AVL/Airfoils/sc20012.dat'

    "AVL HEADER - Plane"
    #Create your plane object. The 'Name' field only displays in AVL, and is not the name of the avl file
    plane = des_plane(Name='TitanWing', Mach=0.8, Sref = 0, bref = 0, cref = 0, Xref=45.478, Yref=0.0, Zref=4.83)
    #The reference parameters are left as zero here      ↑         ↑         ↑
    #This is intentional, as you may choose a surface object to be used to auto-calculate the reference parameters     

    "AVL SURFACE - Surfaces"
    #Create your surfaces.
    wing = des_surf(plane, Name='Wing', YDUPLICATE=0.0, Nspanwise=31)
    hstab = des_surf(plane, Name='Hstab', TRANSLATE=(141.75,0,15), YDUPLICATE=0.0)
    vstab = des_surf(plane, Name='Vstab', TRANSLATE=(131.5,0.0,0.0))
    #Make the plane object use a surface as a reference. Do this to auto-calculate the reference parameters.
    plane.use_reference(wing)
    #AVL documentation suggests that you assign intersecting geometry to use the same COMPONENT number.
    #These parameters may be assigned on object declaration or using .assign()
    #ex: hstab = des_surf(plane, Name='Hstab', TRANSLATE=(141.75,0,15), YDUPLICATE=0.0, COMPONENT=1)
    hstab.assign(COMPONENT=1)
    vstab.assign(COMPONENT=1)

    "AVL SECTION - Sections"
    #Create your sections. Use a G-code-esque naming scheme to allow inserting of sections anywhere
    wingsec000 = des_sec(wing, Xle=0.0, Yle=0.0, Zle=0.0, Chord=cr, Ainc=incidence, AFILE=af1, Nspanwise=21)
    wingsec200 = des_sec(wing, Xle=b/2*np.tan(le_sweep), Yle=b/2, Zle=b/2*np.tan(dihedral), Chord=ct, Ainc=incidence+twist, AFILE=af1) 
    #Another example of resassigning a parameter using .assign()
    wingsec000.assign(Nspanwise=31)

    #horizontal stabilizer sections.
    hstabsec000 = des_sec(hstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=28.8, AFILE=af2)
    hstabsec100 = des_sec(hstab, Xle=31.50934, Yle=45, Zle=0.0, Chord=7.2, AFILE=af2)
    #Use the .interp() function to linearly interpolate between the root and tip
    hstabsec010 = des_sec.interp(s1=hstabsec000, s2=hstabsec100, span=2.5)
    hstabsec090 = des_sec.interp(s1=hstabsec000, s2=hstabsec100, span=40.5)

    #vertical stabilizer sections
    vstabsec000 = des_sec(vstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=35.00161, AFILE=af2)
    vstabsec100 = des_sec(vstab, Xle=27.30110, Yle=0, Zle=42.04, Chord=19.25089, AFILE=af2)

    "AVL CONTROL - Control Surfaces"
    #Keep your control surfaces named "Aileron", "Elevator", and "Rudder" for the best results
    elvstart = des_ctrl(hstabsec010, Cname='Elevator', Xhinge=0.80, SgnDup=1)
    elvend = des_ctrl(hstabsec090, Cname='Elevator', Xhinge=0.80, SgnDup=1)

    "Ordering"
    #You don't need to edit this
    #If using .interp(), this organizes your sections according to their distance from the root. 
    #⚠ AVL breaks if sections are not ordered properly!!
    wing.order(orientation=HORIZONTAL)
    hstab.order(orientation=HORIZONTAL)
    vstab.order(orientation=VERTICAL)

    "Set your filename"
    #Here's where you write avl file. You may decide a custom filename by specifying `filename=(Name)` as a function argument.
    plane.write_to_avl()
