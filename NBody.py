# Sets initial parameters, runs simulations and saves results

from LightSim import *
import time
import numpy as np
from tqdm import tqdm
import math
import sys
import multiprocessing as mp


year = 3600 * 24 * 352.25
AU = 1.49598 * pow(10, 11)
solarMass = 1.988 * pow(10, 30)
parsec = 3.0857 * pow(10, 16)

class Results:
    def __init__(self, bodies, deltaTime, numSteps, runTime):
        self.bodies = bodies
        self.runTime = runTime
        self.time = np.array([i * deltaTime for i in range(numSteps)])



def RunNBody(bodies, deltaTime, numSteps, softeningRadius = 0, saveFactor = None, savePath= ""):
    #returns the results of a single simulation
    ZeroMomentum(bodies)
    CentreOfMass(bodies)

    if "HeavySim" in sys.modules:
        print("initialising vectors...\n")
        InitZeroArrays(bodies, numSteps)

    numCpus = mp.cpu_count()

    print("running simulation...\n")

    time1 = time.time()

    if "HeavySim" in sys.modules:
        for t in tqdm(range(1, numSteps)):
            #print("iteration {} / {}".format(t + 1, numSteps))

            processes = []

            #shared arrays to store data from parallel processing
            positions = [mp.Array("f", 3) for _ in range(len(bodies))]
            velocities = [mp.Array("f", 3) for _ in range(len(bodies))]

            for i in range(numCpus):
                startBodyIndex = round(len(bodies) / numCpus * i)
                endBodyIndex = round(len(bodies) / numCpus * (i+1))

                p = mp.Process(target= RK4, args= [bodies, t, deltaTime, softeningRadius, startBodyIndex, endBodyIndex, positions, velocities])
                p.start()

                processes.append(p)

            for process in processes:
                process.join()

            #save data from shared array
            for i in range(len(bodies)):

                bodies[i].position[t].x = positions[i][0]
                bodies[i].position[t].y = positions[i][1]
                bodies[i].position[t].z = positions[i][2]

                bodies[i].velocity[t].x = velocities[i][0]
                bodies[i].velocity[t].y = velocities[i][1]
                bodies[i].velocity[t].z = velocities[i][2]

    elif "LightSim" in sys.modules:
        if saveFactor != None:
            np.save("{}_0.npy".format(savePath), [[body.position for body in bodies], [body.velocity for body in bodies]],allow_pickle= True)

        for t in tqdm(range(1, numSteps)):

            #print("iteration {} / {}".format(t + 1, numSteps))

            processes = []

            #shared arrays to store data from parallel processing
            positions = [mp.Array("f", 3) for _ in range(len(bodies))]
            velocities = [mp.Array("f", 3) for _ in range(len(bodies))]

            for i in range(numCpus):
                startBodyIndex = round(len(bodies) / numCpus * i)
                endBodyIndex = round(len(bodies) / numCpus * (i+1))

                p = mp.Process(target= RK4, args= [bodies, deltaTime, softeningRadius, startBodyIndex, endBodyIndex, positions, velocities])
                p.start()

                processes.append(p)

            for process in processes:
                process.join()

            #save data from shared array
            for i in range(len(bodies)):
                bodies[i].position.x = positions[i][0]
                bodies[i].position.y = positions[i][1]
                bodies[i].position.z = positions[i][2]

                bodies[i].velocity.x = velocities[i][0]
                bodies[i].velocity.y = velocities[i][1]
                bodies[i].velocity.z = velocities[i][2]
            
            if saveFactor != None:
                if t % saveFactor == 0:
                    np.save("{}_{}.npy".format(savePath, t), [[body.position for body in bodies], [body.velocity for body in bodies]], allow_pickle= True)


    time2 = time.time()
    runTime = time2 - time1

    results = [bodies, [deltaTime, numSteps], runTime, saveFactor]

    return results


def RunMultIts(initialBodies, deltaTime, numSteps, saveFactor= None, savePath = None):
    #returns the results of multiple simulations for varying deltaTime and numSteps
    results = []

    for i in range(len(numSteps)):
        bodies = deepcopy(initialBodies)

        print("\niteration {}\n".format(i +1))

        path = "{}it{}/".format(savePath, i)

        data = RunNBody(bodies, deltaTime[i], numSteps[i], saveFactor= saveFactor, savePath= path)

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
    #sets a power law distributon as the IMF
    np.random.seed(seed= seed)
    
    beta = -alpha +1
    b = beta / (pow(maxMass, beta) - pow(minMass, beta))

    x = np.random.random(len(bodies))

    for i in range(len(bodies)):
        bodies[i].mass = pow((beta * x[i]) / b + pow(minMass, beta), 1/beta)


    totalBodiesMass = sum([body.mass for body in bodies])

    if clusterMass != None:
        for body in bodies:
            body.mass *= clusterMass / totalBodiesMass


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
    r = radius * np.power(x3, 1/3)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    for i in range(len(bodies)):
        bodies[i].position = Vector3(x[i], y[i], z[i])


def SetVirialVelocityDist(bodies, seed, clusterMass, clusterRadius, velocityModifier = 1):
    #for uniform position distribution
    density = clusterMass / (4/3 * np.pi * pow(clusterRadius, 3))

    #random direction of velocity along sphere of radius = body.position.Magnitude()

    np.random.seed(seed= seed + 3)
    x1 = np.random.random(len(bodies))

    np.random.seed(seed= seed + 4)
    x2 = np.random.random(len(bodies))

    theta = np.arccos(1 - 2 * x1)
    phi = 2 * np.pi * x2

    enclosedMass = np.array([4 / 3 * np.pi * pow(bodies[i].position.Magnitude(), 3) * density for i in range(len(bodies))])
    velocityMag = np.array([pow(G * enclosedMass[i] / bodies[i].position.Magnitude(), 1/2) * velocityModifier for i in range(len(bodies))])

    x = velocityMag * np.sin(theta) * np.cos(phi)
    y = velocityMag * np.sin(theta) * np.sin(phi)
    z = velocityMag * np.cos(theta)

    for i in range(len(bodies)):
        bodies[i].velocity = Vector3(x[i], y[i], z[i])


if __name__ == "__main__":
    endTime = year * 1 * pow(10, 6)
    """
    # Use below for single run
    """
    deltaTime = year * 1 * pow(10, 0)
    numSteps = round(endTime / deltaTime)
    """
    # Use below for multiple runs
    # """
    # startDeltaTime = year * 1 * pow(10, 1)
    # endDeltaTime = year * 1 * pow(10, 3)
    # numIterations = 3
    # deltaTime, numSteps = GetLogTimeArrays(endTime, numIterations, startDeltaTime, endDeltaTime)
    """
    """

    numBodies = 1 * pow(10, 1)
    bodies = [Body.Zero() for _ in range(numBodies)]

    seed = 3

    clusterAlpha = 2.35
    clusterMass = 1.3 * pow(10, 6) * solarMass
    clusterRadius = 6 * parsec

    minMass = 0.1 * solarMass
    maxMass = 100 * solarMass

    velocityModifier = 1

    softeningRadius = 0

    SetIMFPowerLaw(bodies, seed, minMass, maxMass, clusterAlpha, clusterMass)
    SetUniformPositionDist(bodies, seed, clusterRadius)
    SetVirialVelocityDist(bodies, seed, clusterMass, clusterRadius, velocityModifier= velocityModifier)

    savePath = "ComGPE_dtE0_Seed{}/".format(seed)
    saveFactor = 100

    print(savePath)

    if type(deltaTime) == list:
        results = RunMultIts(bodies, deltaTime, numSteps, savePath= savePath, saveFactor= saveFactor)
    
    elif type(deltaTime) == float:
        results = RunNBody(bodies, deltaTime, numSteps, saveFactor= saveFactor, savePath= savePath, softeningRadius= softeningRadius)

    print("\nsaving results...\n")

    if "LightSim" in sys.modules:
        np.save("{}_end.npy".format(savePath), results, allow_pickle= True)
    
    elif "HeavySim" in sys.modules:
        np.save("{}.npy".format(savePath), results, allow_pickle= True)
    
    print("save complete\n")




