from computation import nozzle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Inputs
Me = 2.286

# Constants
GAMMA = 1.4
n_list = [8,16,32,64,128]
R = .5
TOL = 1e-12

def polyfitfunction(x,a,b,c,d):
    return a*np.exp(-b*x+c)+d

if __name__ == '__main__':

    area_ratio_list = np.zeros(len(n_list))

    for i, n in enumerate(n_list):
        area_ratio = nozzle(GAMMA, Me, n, R, i)
        area_ratio_list[i] = (area_ratio)

    plt.figure(i+2)

    exact_AR = (2.4/2) ** (-3) * ((1+.2*Me**2) ** (3))/Me
    error = np.abs(exact_AR - area_ratio_list)/exact_AR * 100

    plt.scatter(n_list, error, s=30, color='r', zorder=15, label="Data")

    popt, _ = curve_fit(polyfitfunction, np.array(n_list),error, bounds=(0,2))
    a,b,c,d = popt

    x_error = np.linspace(min(n_list), max(n_list), 1000)
    y_error = polyfitfunction(x_error,a,b,c,d)

    plt.plot(x_error,y_error, label="Exponential Interpolation Fit")
    plt.xlabel('Number of Initial Points')
    # Plots y label
    plt.ylabel('True Percentage Relative Error [%]')
    # Plots title
    plt.title('True Percentage Relative Error vs. Number of Initial Points')
    plt.legend(loc="upper right")

    plt.show()
