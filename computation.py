import numpy as np
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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

def solve_M(z, gamma, Me):
    func = lambda M: np.sqrt((gamma+1)/(gamma-1)) *\
                     np.degrees(np.arctan(np.sqrt((gamma-1)*(M ** 2 -1)/(gamma+1))))\
                     - np.degrees(np.arctan(np.sqrt((M ** 2 -1)))) - z
    tau_initial_guess = Me
    return fsolve(func, tau_initial_guess)

def sline(L,xL,yL,R,xR,yR):

    x = (yR-yL+L*xL-R*xR)/(L-R)
    y = yR+R*(x-xR)
    return x,y

def nozzle(gamma, Me, n, R):

    R = np.zeros([(n+2), n])
    L = np.zeros([(n+2), n])

    d_max = nu(Me, gamma)/2
    d_inc = d_max/n

    d_init = np.linspace(d_inc,d_max, n)

    v = d_init.copy()
    L[0,:] = d_init - v
    R[0,:] = d_init + v

    # Field Points

    # Right Points

    for i in range(n):
        R[1:(i+2),i:(i+1)] = np.ones([(i+1),1]) * R[0,i]

    # Wall Points
    R[(R==0)] = R[0,-1]

    # Left

    for i in range(n):
        L[(i+1):(i+2), i:(n)] = np.ones([1, (n-i)]) * R[0, i]

    # Wall Points
    for i in range(n):
        L[(i+2):(i+3), 0:(i+1)] = np.ones([1, (i+1)]) * R[0, i]

    d = (R-L)/2
    v = (R+L)/2

    vector_solve = np.vectorize(solve_M)

    M = vector_solve(v, gamma, Me)

    u = np.degrees(np.arcsin(1/M))

    r = 0.5

    R_slope = np.zeros([d.shape[0]-1, d.shape[1]])

    R_slope[0,:] = np.tan(np.radians((d[0,:]-u[0,:])))

    R_slope[1:,:] = np.tan(np.radians(d[1:-1,:] - u[1:-1,:] + d[:-2,:] - u[:-2,:])/2)

    R_slope_diag_msk = np.eye(R_slope.shape[0], R_slope.shape[1], k=-2, dtype=bool)
    R_slope[R_slope_diag_msk] = np.tan(np.radians(d[:-1,:][R_slope_diag_msk] - u[:-1,:][R_slope_diag_msk]))

    R_slope_diag_msk = np.eye(R_slope.shape[0], R_slope.shape[1], k=-1, dtype=bool)
    R_slope[R_slope_diag_msk] = 0

    L_slope = np.zeros(R.shape)

    L_slope[1:,1:] = np.tan(np.radians((d[1:,1:] + u[1:,1:] + d[1:,:-1] + u[1:,:-1])/2))
    L_slope[1:,0] = np.tan(np.radians((d[1:,0] + u[1:,0] + d[:-1, -1] + u[:-1,-1])/2))

    L_slope_diag_msk = np.eye(L_slope.shape[0], L_slope.shape[1], k=-2, dtype=bool)
    L_slope[L_slope_diag_msk] = np.tan(np.radians(d_max-d_init))

    x = np.zeros(R.shape)
    y = np.zeros(R.shape)

    x[0,:] = r * np.sin(np.radians(d_init))
    y[0,:] = 2 * r - r * np.cos(np.radians(d_init))

    # x[1:,:] = sline()
    # y[1:,:] = sline()

    # Can't use vectorization since some cells rely on previous cells

    for i in range(1,n+2):
        for j in range(n):
            if i == 2 and j == 0:
                x[i,j],y[i,j] = sline(L_slope[i-1, -1], x[i-1,-1], y[i-1,-1], np.tan(np.radians(d_max-d_inc)), x[0,-1], y[0,-1])
            elif (j)== (i-2):
                x[i,j], y[i,j] = sline(L_slope[i, j-1], x[i,j-1], y[i,j-1], L_slope[i-1,j-1],x[i-1,j-1],y[i-1,j-1])
            elif not (j) == (i-1):
                if j == 0:
                    x[i,j], y[i,j] = sline(L_slope[i-1,-1],x[i-1,-1],y[i-1,-1],R_slope[i-1,0],x[i-1,0],y[i-1,0])
                else:
                    x[i,j], y[i,j] = sline(L_slope[i,j-1],x[i,j-1],y[i,j-1],R_slope[i-1,j],x[i-1,j],y[i-1,j])
            else:
                x[i,j], y[i,j] = sline(0,0,0, R_slope[i-1,j], x[i-1,j], y[i-1,j])

    fig, axs = plt.subplots(2)

    lw = .85
    rgb = [.85, .85, .85]

    for i in range(n):
        axs[1] = plt.plot(x[0:(i+2),i], y[0:(i+2),i], color = rgb, linewidth = lw)
        axs[1] = plt.plot(x[2:,i], y[2:,i], color = rgb, linewidth = lw)
        axs[1] = plt.plot(np.concatenate((x[i+1, i:], x[i+2,:i+1])), np.concatenate((y[i+1, i:], y[i+2,:i+1])), color = rgb, linewidth = lw)

    axs[1] = plt.plot(np.concatenate((x[0, :], x[L_slope_diag_msk])), np.concatenate((y[0, :], y[L_slope_diag_msk])), color = [0,0,0], linewidth = (lw))
    axs[1] = plt.scatter(x, y, c=M, cmap='viridis', s = 15, zorder = 15)
    axs[1] = plt.colorbar(orientation="vertical")
    axs[1].ax.set_title('M')


    axs[1] = plt.grid()
    axs[1] = plt.xlabel('x/L')
    axs[1] = plt.ylabel('y/L')
    axs[1] = plt.title('M = ' + str(Me) + ", N = " + str(n) + ", Area Ratio = " + "{:.3f}".format(2*y[-1,-1]))
    plt.show()

    l = 0

    return

















