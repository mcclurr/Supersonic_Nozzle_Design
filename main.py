
from computation import nozzle


# Inputs
Me = 2.286

# Constants
GAMMA = 1.4
n = [8, 16]
R = .5

if __name__ == '__main__':

    area_ratio = []

    for i in n:
        area_ratio.append(nozzle(GAMMA, Me, n, R))




