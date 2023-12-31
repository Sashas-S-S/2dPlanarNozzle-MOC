# 2dPlanarNozzle-MOC
This is a matlab Code for generating the contour of a Supersonic 2d Planar Nozzle using Method of Characteristics based on outlet Mach number, adiabatic index of the operational fluid, number of Characteristic Lines and throat diameter. The wall contour generated is only for the divergent section of the nozzle assuming a centered Prandtl Meyer expansion fan at the nozzle throat. This assumption generates a minimum length divergent section contour needed for the nozzle. Based on other parameters like chamber pressure, temperature and ambient pressure, temperature other characteristics of the nozzle like thrust and expansion ratio are also calculated using analytical relations. 

# Valid for Ma < 9

# Materials Used
All the assumptions for the calculations and material used for the code are below:
1. Modern Compressible Flow, with Historical Perspective,Anderson, J.D, 9780070016545,McGraw-Hill series in aeronautical and aerospace engineering.
2. NPTEL, Advanced Gas Dynamics by Dr.Rinku Mukherjee, IITM.
3. NPTEL, Aerospace Propulsion by Dr. P.A. Ramakrishna, IITM.

# Run code from methodOfCharateristics.m
Export coordinates of the nozzle from the variables xcoords,ycoords and use them to import into any CAD software.
