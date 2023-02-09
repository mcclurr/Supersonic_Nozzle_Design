
from computation import nozzle
import matplotlib.pyplot as plt


# Inputs
Me = 2.286

# Constants
GAMMA = 1.4
n_list = [8,16,32,64,128]
R = .5

if __name__ == '__main__':

    for i, n in enumerate(n_list):
        fig = nozzle(GAMMA, Me, n, R, i)

    plt.show()
