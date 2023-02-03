import numpy as np
import math

###################################### Inputs ######################################

# gamma: Heat Coefficient Ratio
# Me: Exit Mach Number
# N: number of initial points
# R: throat expansion curvature radius with respect to throat height

###################################### Outputs #####################################

# M: Mach Number Distribution

######################### Variables Used for Generating the Shape ##################

# delta: d
# Prandtl-Meyer Angle: v
# C- Characteristic Line: L
# C+ Characteristic Line: R
# Mach: M
# Mach Angle: u
# Slope of L Lines: Lslope
# Slope of R Lines: Rslope
# x
# y

##################################### Functions ###################################

def nu(M, gamma):

    theta = math.sqrt((gamma+1)/(gamma-1)) \
        * math.degrees(math.atan(math.sqrt((gamma-1)*(M ** 2-1)/(gamma+1))))\
        - math.degrees(math.atan(math.sqrt(M ** 2 - 1)))
    return theta



def nozzle(gamma, Me, n, R):

    d_max = nu(Me, gamma)/2
    d_inc = d_max/n

    d = np.linspace(d_inc,d_max, n)

    v = d.copy()
    L = d - v
    R = d + v
    


    l = 0

    return x

















