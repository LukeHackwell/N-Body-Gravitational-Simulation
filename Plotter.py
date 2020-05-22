# Plot results from NBody.py

from matplotlib import pyplot as plt
import numpy as np
from HeavySim import *
import LightSim
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

year = 3600 * 24 * 352.25
AU = 149598 * pow(10, 6)
solarMass = 1.988 * pow(10, 30)
G = 6.6743 * pow(10, -11)
parsec = 3.0857 * pow(10, 16)



def FilterByGPE(bodies, filterFactor, index):
    gPE = CalcGPE(bodies)
    gPE = [gPE[i][-1] for i in range(len(bodies))]
    minGPE = min(gPE)

    filteredBodies = []
    for i in range(len(gPE)):
        if gPE[i] < (minGPE / filterFactor):
            filteredBodies.append(bodies[i])

    return filteredBodies



def GetRadiiTimeMassFromMultSave(savePath, endSavePath= None, saveIndex= None):
    if endSavePath != None:
        endRes = np.load(endSavePath, allow_pickle= True)
    else:
        endRes = np.load("{}_end.npy".format(savePath), allow_pickle= True)

    if endSavePath == None:
        deltaTime = endRes[1][0]
        numSteps = endRes[1][1]
        saveFactor = endRes[3]


        lenBodies = len(endRes[0])

        time = [0]
        for t in range(1, numSteps):
            if t % saveFactor == 0:
                time.append(t)

        mass = [0 for _ in range(lenBodies)]
        for i in range(len(endRes[0])):
            print(endRes[0][i])
            mass[i] = endRes[0][i].mass

        radii = [[] for _ in range(lenBodies)]
        for t in time:
            data = np.load("{}_{}.npy".format(savePath, t), allow_pickle= True)
            for i in range(lenBodies):
                r = data[0][i].Magnitude()
                radii[i].append(r)

        return radii, time, mass

    if type(endRes[0][0]) == list:
        deltaTime = endRes[saveIndex][1][0]
        numSteps = endRes[saveIndex][1][1]
        saveFactor = 100

        lenBodies = len(endRes[saveIndex][0])

        time = [0]

        for t in range(1, numSteps):
            if t % saveFactor == 0:
                time.append(t)

        mass = [0 for _ in range(lenBodies)]
        for i in range(len(endRes[saveIndex][0])):
            mass[i] = endRes[saveIndex][0][i].mass

        radii = [[] for _ in range(lenBodies)]

        for t in time:
            data = np.load("{}_{}.npy".format(savePath, t), allow_pickle= True)
            for i in range(lenBodies):
                r = data[0][i].Magnitude()
                radii[i].append(r)

        return radii, time, mass


def GetPositionsFromMultSave(savePath, endSavePath= None, saveIndex= None):
    if endSavePath != None:
        endRes = np.load(endSavePath, allow_pickle= True)
    else:
        endRes = np.load("{}_end.npy".format(savePath), allow_pickle= True)

    if endSavePath == None:
        deltaTime = endRes[1][0]
        numSteps = endRes[1][1]
        saveFactor = endRes[3]


        lenBodies = len(endRes[0])

        time = [0]

        for t in range(1, numSteps):
            if t % saveFactor == 0:
                time.append(t)

        mass = [0 for _ in range(lenBodies)]
        for i in range(len(endRes[0])):
            mass[i] = endRes[0][i].mass

        positions = [[] for _ in range(lenBodies)]

        for t in time:
            data = np.load("{}_{}.npy".format(savePath, t), allow_pickle= True)
            for i in range(lenBodies):
                position = data[0][i]
                positions[i].append(position)

        return positions, time, mass

    if type(endRes[0][0]) == list:
        deltaTime = endRes[saveIndex][1][0]
        numSteps = endRes[saveIndex][1][1]
        saveFactor = 100

        lenBodies = len(endRes[saveIndex][0])

        time = [0]

        for t in range(1, numSteps):
            if t % saveFactor == 0:
                time.append(t) / pow(10, 6) / year

        mass = [0 for _ in range(lenBodies)]
        for i in range(len(endRes[saveIndex][0])):
            mass[i] = endRes[saveIndex][0][i].mass

        positions = [[] for _ in range(lenBodies)]

        for t in time:
            data = np.load("{}_{}.npy".format(savePath, t), allow_pickle= True)
            for i in range(lenBodies):
                position = data[0][i]
                positions[i].append(position)

        return positions, time, mass


def GetBodiesTimeFromMultSave(savePath, endSavePath= None, saveIndex= None):
    positions, time, mass = GetPositionsFromMultSave(savePath, endSavePath= endSavePath, saveIndex= saveIndex)

    bodies = [Body.Zero() for _ in range(len(positions))]
    InitZeroArrays(bodies, len(positions[0]))

    for i in range(len(bodies)):
        for t in range(len(positions[i])):
            bodies[i].position[t] = positions[i][t]

        bodies[i].mass = mass[i]

    return bodies, time




"""
Get results that compare final radius and dt
"""
# deltaTime = []
# numSteps = []
# endTime = []

# radii = [[] for _ in range(len(results[0][0]))]

# i = 0
# for iteration in results:
#     deltaTime.append(iteration[1][0] / year)
#     numSteps.append(iteration[1][1])
#     endTime.append(iteration[1][0] * iteration[1][1])


# for i in range(len(results)):
#     for j in range(len(results[i][0])):
#         radii[j].append(results[i][0][j].position.Magnitude() / parsec)

"""
 plot final radius for each iteration
"""
# for i in range(len(radii)):
#     plt.plot(deltaTime, radii[i], label= i)

# plt.xlabel("Delta Time (years)")
# plt.ylabel("Final Radius (pc)")
# plt.title("Final Radius of Bodies vs Delta Time")
# plt.xscale("log")
# plt.legend()
# plt.show()

"""
 find percentage difference in all bodies
"""

# unboundIndicies = []
# boundRadii = []

# for i in range(len(radii)):
#     if i not in unboundIndicies:
#         boundRadii.append(radii[i])

# radii = boundRadii

# print("delta times: {}".format(deltaTime))

# index1 = 0
# index2 = 2
# print("min delta time: {}, max delta time: {} \n".format(deltaTime[index1], deltaTime[index2]))

# percentage = [0 for _ in range(len(radii))]
# for i in range(len(radii)):
#     percentage[i] = (radii[i][index2] - radii[i][index1]) / radii[i][index1]

# sumPercPerBody = sum(percentage) / len(radii) * 100

# print("percentage change per body: {}\n".format(sumPercPerBody))


"""
 print max percentage difference in final radii
""""
# maxFractionalDifference = 0
# maxDeltaTimeIndex = 0
#
# maxDeltaTimeIndex = len(deltaTime) - 0
#
# for j in range(len(radii)):
#
#     minRadii = radii[j][0]
#     maxRadii = radii[j][0]
#
#     for i in range(maxDeltaTimeIndex):
#         if radii[j][i] < minRadii:
#             minRadii = radii[j][i]
#
#         if radii[j][i] > maxRadii:
#             maxRadii = radii[j][i]
#
#     fractionalDifference = (maxRadii - minRadii) / radii[j][0]
#
#     if fractionalDifference > maxFractionalDifference:
#         maxFractionalDifference = fractionalDifference
#
# print("max percentage difference = {}".format(maxFractionalDifference * 100))
# print("corresponding max delta time = {} years".format(deltaTime[-1]))




"""
Plotting info at each time step:
"""
# bodyIndex = 0
"""
 plot 3D track of one body's positions
"""
# fig = plt.figure()
# ax = fig.add_subplot(111, projection= "3d")
# Axes3D.plot(ax, positions[bodyIndex][0], positions[bodyIndex][1], zs =positions[bodyIndex][2])
# plt.show()

"""
 plot distance from CM for one body
"""
#distance = np.array([bodies[bodyIndex].position[i].Magnitude() for i in range(len(bodies[bodyIndex].position))]) / parsec

#plt.plot(time, distance)
#plt.xlabel("time (years)")
#plt.ylabel("distance from CM (parsec)")
#plt.show()

"""
 plot 3d track of all body positions
"""
# positions, velocities = VerticalArray(bodies, parsec)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection= "3d")
# for i in range(len(positions)):
#     Axes3D.plot(ax, positions[i][0], positions[i][1], zs =positions[i][2], label= str(i))

# plt.title("Positions of bodies (parsec). Legend: Body Index")
# plt.legend()
# plt.show()

"""
 plot distance from centre for all bodies
"""
# distance = [0 for i in range(len(bodies))]
# for i in range(len(bodies)):
#     distance[i] = np.array([bodies[i].position[j].Magnitude() for j in range(len(bodies[i].position))]) / parsec

# for i in range(len(distance)):
#     plt.plot(time, distance[i], label= str(i))

# plt.xlabel("time (yr)")
# plt.ylabel("distance from CM (parsec)")
# plt.legend()
# plt.show()


"""
 plot GPE vs time for all bodies
"""
# gPE = CalcGPE(bodies)

# for i in range(len(gPE)):
#     plt.plot(time, gPE[i], label= str(i))

# plt.xlabel("time (yr)")
# plt.ylabel("GPE (J)")
# plt.legend()
# plt.show()

"""
 plot GP vs time for all bodies
"""
# gP = CalcGPE(bodies)

# for i in range(len(gP)):
#     plt.plot(time, gP[i], label= str(i))

# plt.xlabel("time (yr)")
# plt.ylabel("GP (J/kg)")
# plt.legend()
# plt.show()

