

from computation import nozzle
import matplotlib.pyplot as plt
import numpy as np


# Inputs
Me = 2.286

# Constants
GAMMA = 1.4
n_list = [8,16,32,64,128]
R = .5

if __name__ == '__main__':

    area_ratio_list = np.zeros(len(n_list))

    for i, n in enumerate(n_list):
        area_ratio = nozzle(GAMMA, Me, n, R, i)
        area_ratio_list[i] = (area_ratio)

    plt.figure(i+2)

    exact_AR = (2.4/2) ** (-3) * ((1+.2*Me**2) ** (3))/Me
    error = np.abs(exact_AR - area_ratio_list)/exact_AR * 100

    plt.scatter(n_list, error, s=20, facecolors='none', edgecolors='g')

    l = 0


    plt.show()
