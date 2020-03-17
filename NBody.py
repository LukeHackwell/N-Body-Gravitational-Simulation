from HeavySim import *
import time
import numpy as np
from tqdm import tqdm
import math
import sys


year = 3600 * 24 * 352.25
AU = 1.49598 * pow(10, 11)
solarMass = 1.988 * pow(10, 30)
parsec = 3.0857 * pow(10, 16)

def RunNBody(bodies, deltaTime, numSteps, function, softeningRadius = 0):
    #returns the results of a single simulation

    ZeroMomentum(bodies)
    CentreOfMass(bodies)
    InitSmartAcceleration(bodies)

    if "HeavySim" in sys.modules:
        InitZeroArrays(bodies, numSteps)

    time1 = time.time()

    if "HeavySim" in sys.modules:
        for t in tqdm(range(1, numSteps)):
            function(bodies, t, deltaTime, softeningRadius)

    else:
        for _ in tqdm(range(1, numSteps)):
            function(bodies, deltaTime, softeningRadius)

    time2 = time.time()
    runTime = time2 - time1

    results = [bodies, [deltaTime, numSteps], runTime]

    return results


def RunMultIts(initialBodies, deltaTime, numSteps, function):
    #returns the results of multiple simulations for varying deltaTime and numSteps

    results = []

    for i in range(len(numSteps)):
        bodies = deepcopy(initialBodies)

        data = RunNBody(bodies, deltaTime[i], numSteps[i], function)

        results.append(data)

    return results


def GetLinearTimeArrays(endTime, numIterations, startDeltaTime, endDeltaTime):
    # returns a linear range of delta time and num steps

    deltaTime = [0 for i in range(numIterations)]
    numStepsArray = [0 for i in range(numIterations)]
    deltaTimeOffset = (endDeltaTime - startDeltaTime) / numIterations

    for i in range(numIterations):
        deltaTime[i] = startDeltaTime + i * deltaTimeOffset
        numStepsArray[i] = round(endTime / deltaTime[i])

    return deltaTime, numStepsArray


def GetLogTimeArrays(endTime, numIterations, startDeltaTime, endDeltaTime):
    #returns a logarithmic range of delta time and num steps

    deltaTime = [0 for i in range(numIterations)]
    numStepsArray = [0 for i in range(numIterations)]
    N = math.log10(endDeltaTime / startDeltaTime)

    for i in range(numIterations):
        deltaTime[i] = startDeltaTime * pow(10, N * i / numIterations)
        numStepsArray[i] = round(endTime / deltaTime[i])

    deltaTime.append(startDeltaTime * pow(10, (N * numIterations) / numIterations))
    numStepsArray.append(round(endTime / deltaTime[-1]))

    return deltaTime, numStepsArray


def SetIMFPowerLaw(bodies, seed, minMass, maxMass, alpha, clusterMass= None):
    #sets a power law distributon for the masses

    np.random.seed(seed= seed)
    
    beta = -alpha +1
    b = beta / (pow(maxMass, beta) - pow(minMass, beta))

    x = np.random.random(len(bodies))

    for i in range(len(bodies)):
        bodies[i].mass = pow((beta * x[i]) / b + pow(minMass, beta), 1/beta)

    
    if clusterMass != None:
        totalBodiesMass = sum([body.mass for body in bodies])
        
        for body in bodies:
            body.mass *= clusterMass / totalBodiesMass / len(bodies)


def SetUniformPositionDist(bodies, seed, radius):
    #Sets bodies positions uniformly within cluster radius
    #x is a np random number array. They have different seeds so that phi, theta, and r have different random values

    np.random.seed(seed= seed)
    x1 = np.random.random(len(bodies))

    np.random.seed(seed= seed+1)
    x2 = np.random.random(len(bodies))

    np.random.seed(seed= seed+2)
    x3 = np.random.random(len(bodies))

    theta = np.arccos(1 - 2 * x1)
    phi = 2 * np.pi * x2
    r = radius * np.power(4/3 * np.pi * x3, 1/3)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    for i in range(len(bodies)):
        bodies[i].position = Vector3(x[i], y[i], z[i])


def SetVelocityDistribution(bodies, seed):

    pass
        
        
        
        
        
        
        
    


if __name__ == "__main__":
    endTime = year * 1

    """
    Use below for single run
    """
    deltaTime = 3600
    numSteps = round(endTime / deltaTime)

    """
    Use below for multiple runs
    """
    # startDeltaTime = 3600
    # endDeltaTime = 3600 * 1000
    # numIterations = 5
    # deltaTime, numSteps = GetLogTimeArrays(endTime, numIterations, startDeltaTime, endDeltaTime)



    numBodies = pow(10, 1)
    bodies = [Body.Zero() for _ in range(numBodies)]

    seed = 1

    clusterAlpha = 2.49
    clusterMass = 2.3 * pow(10, 6) * solarMass
    clusterRadius = parsec


    SetIMFPowerLaw(bodies, seed, 0.1 * solarMass, 10 * solarMass, clusterAlpha, clusterMass= clusterMass)
    SetUniformPositionDist(bodies, seed, clusterRadius)


    print("running simulation...")
    results = RunNBody(bodies, deltaTime, numSteps, RK4)
    print("saving results...")
    np.save("Results.npy", results, allow_pickle= True)







