from matplotlib import pyplot as plt
import numpy as np
from HeavySim import *
from mpl_toolkits.mplot3d import Axes3D
import sys

year = 3600 * 24 * 352.25
AU = 149598 * pow(10, 6)
solarMass = 1.988 * pow(10, 30)
G = 6.6743 * pow(10, -11)
parsec = 3.0857 * pow(10, 16)






results = np.load("Results.npy", allow_pickle= True)

bodies = results[0]

positions, velocities = VerticalArray(bodies, parsec)

deltaTime = results[1][0]
numSteps = results[1][1]
endTime = deltaTime * numSteps

time = np.array([i * deltaTime for i in range(numSteps)]) / year


"""
LightSim
"""
if "LightSim" in sys.modules:
    pass

    # #Plot final positions of all bodies
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection= "3d")
    # Axes3D.scatter3D(ax, positions[0], positions[1], zs =positions[2])
    # plt.show()




"""
HeavySim
"""
if "HeavySim" in sys.modules:
    bodyIndex = 0

    # # plot 3D track of one body's positions
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection= "3d")
    # Axes3D.scatter3D(ax, positions[bodyIndex][0], positions[bodyIndex][1], zs =positions[bodyIndex][2])
    # plt.show()




    # # plot distance from CM
    # distance = np.array([bodies[bodyIndex].position[i].Magnitude() for i in range(len(bodies[bodyIndex].position))]) / parsec
    #
    # plt.plot(time, distance)
    # plt.xlabel("time (years)")
    # plt.ylabel("distance from CM (parsec)")
    # plt.show()








