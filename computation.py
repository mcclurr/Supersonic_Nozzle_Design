import numpy as np
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy

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

    tau_initial_guess = 1
    return fsolve(func, tau_initial_guess)

def sline(L,xL,yL,R,xR,yR):

    x = (yR-yL+L*xL-R*xR)/(L-R)
    y = yR+R*(x-xR)
    return x,y

def plot_figs(x,y,M,n,figNum,r,Me):

    diag_msk = np.eye(x.shape[0], x.shape[1], k=-2, dtype=bool)

    # Initializes Figure
    plt.figure(figNum+1)
    # Initializes top Subplot (plots data points)
    plt.subplot(211)

    # line width
    lw = .85
    # Plots lines as Grey
    rgb = [.85, .85, .85]

    # Plots lines between points
    for i in range(n):
        plt.plot(x[0:(i+2),i], y[0:(i+2),i], color = rgb, linewidth = lw)
        plt.plot(x[2:,i], y[2:,i], color = rgb, linewidth = lw)
        plt.plot(np.concatenate((x[i+1, i:], x[i+2,:i+1])), np.concatenate((y[i+1, i:], y[i+2,:i+1])), color = rgb, linewidth = lw)

    # Plots black line to outline edge of nozzle
    plt.plot(np.concatenate((x[0, :], x[diag_msk])), np.concatenate((y[0, :], y[diag_msk])), color = [0,0,0], linewidth = (lw), zorder = 25)
    # Plots scatter plot of x and y points colored by Mach Value
    plt.scatter(x, y, c=M, cmap='viridis', s = 10, zorder = 15)
    # Plots colorbar of Mach
    cbar = plt.colorbar(orientation="vertical")
    # Sets title of 'Mach'
    cbar.ax.set_title('Mach')

    # Set to makes axis equal
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_ylim(top=1.5)
    ax.set_xlim(0,x.max())

    # Plots grid
    plt.grid()
    # Plots x label
    plt.xlabel('x/L')
    # Plots y label
    plt.ylabel('y/L',rotation=0)
    # Plots title
    plt.title('Mach = ' + str(Me) + ", N = " + str(n) + ", Area Ratio = " + "{:.3f}".format(2*y[-1,-1]))


    # Starts computation to plot second subplot

    # Creates X,Y meshgrid
    X,Y = np.meshgrid(np.linspace(0,x.max(),1000), np.linspace(0,y.max(),1000))

    # pulls x,y,z values to be interpolated
    x_intrp = np.concatenate((np.append([0,0],x[0,:]), np.append(np.asarray(x).ravel(), x.max())))
    y_intrp = np.concatenate((np.append([0,r],y[0,:]), np.append(np.asarray(y).ravel(), 0)))
    z_intrp = np.concatenate((np.append([1,1],M[0,:]), np.append(np.asarray(M).ravel(), Me)))

    # Creates Z meshgrid values
    xi = (x_intrp, y_intrp)
    xx = (np.asarray(X), np.asarray(Y))
    Z = scipy.interpolate.griddata(xi, z_intrp, xx)

    # Starts plotting second subplot
    plt.subplot(212)
    plt.pcolormesh(X,Y,Z, cmap='viridis')

    plt.plot(np.concatenate((x[0, :], x[diag_msk])), np.concatenate((y[0, :], y[diag_msk])), color=[0, 0, 0], linewidth=(lw), zorder = 25)

    # Plots colorbar
    cbar = plt.colorbar(orientation="vertical")
    cbar.ax.set_title('Mach')

    # Sets Axis Equal
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_ylim(top=1.5)
    ax.set_xlim(0,x.max())

    # Plots Grid
    plt.grid()
    plt.xlabel('x/L')
    plt.ylabel('y/L', rotation=0)
    plt.title('Mach = ' + str(Me) + ", N = " + str(n) + ", Area Ratio = " + "{:.3f}".format(2*y[-1,-1]))

    plt.figure(figNum+1).tight_layout()

    return

def nozzle(gamma, Me, n, r, figNum):

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

    # Calls on fsolve to numerically solve for the Mach
    vector_solve = np.vectorize(solve_M)
    M = vector_solve(v, gamma, Me)

    u = np.degrees(np.arcsin(1/M))

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

    plot_figs(x,y,M,n,figNum,r,Me)

    area_ratio = 2*y[-1,-1]

    return area_ratio
