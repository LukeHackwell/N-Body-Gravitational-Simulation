from copy import deepcopy
import numpy as np

G = 6.6743 * pow(10, -11)


def Zeros(n):
    #creates an array of zeros of length n
    return [0 for _ in range(n)]


class Vector3:
    #stores 3D vectors


    def __init__(self, x, y, z):
        self.x = x
        self.y = y        
        self.z = z

    def __add__(self, vector2):
        #returns self + vector2
        return Vector3(self.x + vector2.x, self.y + vector2.y, self.z + vector2.z)
    
    def __sub__(self, vector2):
        #returns subtracts vector2 from self
        return Vector3(self.x - vector2.x, self.y - vector2.y, self.z - vector2.z)

    def __mul__(self, factor):
        #returns self multiplied by scalar factor
        return Vector3(self.x * factor, self.y * factor, self.z * factor)

    def __truediv__(self, factor):
        #returns self divided by scalar factor
        return Vector3(self.x / factor, self.y / factor, self.z / factor)

    def __iadd__(self, vector2):
        return Vector3(self.x + vector2.x, self.y + vector2.y, self.z + vector2.z)

    def __isub__(self, vector2):
        return Vector3(self.x - vector2.x, self.y - vector2.y, self.z - vector2.z)

    def __imul__(self, factor):
        return Vector3(self.x * factor, self.y * factor, self.z * factor)

    def __itruediv__(self, factor):
        return Vector3(self.x / factor, self.y / factor, self.z / factor)

    def __str__(self):
        return "[{}, {}, {}]".format(self.x, self.y, self.z)


    def Magnitude(self):
        #return the magnitude of self
        return pow(pow(self.x, 2) + pow(self.y, 2) + pow(self.z, 2), 0.5)


    @staticmethod
    def Zero():
        #return a zero vector
        return Vector3(0, 0, 0)



class Vector3Array:
    #stores a list of 3D vectors. Used to make arithmetic operations easy


    def __init__(self, array):
        self.array = array

    def __add__(self, other):
        return Vector3Array([self.array[i] + other.array[i] for i in range(len(self.array))])

    def __sub__(self, other):
        return Vector3Array([self.array[i] - other.array[i] for i in range(len(self.array))])

    def __mul__(self, factor):
        return Vector3Array([self.array[i] * factor for i in range(len(self.array))])

    def __truediv__(self, factor):
        return Vector3Array([self.array[i] / factor for i in range(len(self.array))])

    def __iadd__(self, other):
        return Vector3Array([self.array[i] + other.array[i] for i in range(len(self.array))])

    def __isub__(self, other):
        return Vector3Array([self.array[i] - other.array[i] for i in range(len(self.array))])

    def __imul__(self, factor):
        return Vector3Array([self.array[i] * factor for i in range(len(self.array))])

    def __itruediv__(self, factor):
        return Vector3Array([self.array[i] / factor for i in range(len(self.array))])

    def __getitem__(self, index):
        return self.array[index]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __str__(self):
        string = ""
        for vector in self.array:
            string += str(vector) + "\n"
        return string

    def Sum(self):
        sum = Vector3.Zero()
        for element in self.array:
            sum += element

        return sum

    
class Body:
    #stores information about a simulated body


    def __init__(self, position, velocity, mass):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.totalAcceleration = 0

    @staticmethod
    def Zero():
        return Body(Vector3(0, 0, 0), Vector3(0, 0, 0), 0)

    def ZerosPosVel(self, n):
        #make a horizontal zeros array for position and velocity
        initialPosition = self.position
        initialVelocity = self.velocity

        zeros = [Vector3.Zero() for _ in range(n)]

        self.position = deepcopy(zeros)
        self.velocity = deepcopy(zeros)

        self.position[0] = initialPosition
        self.velocity[0] = initialVelocity


def InverseSq(r1, r2, softeningRadius):
    #returns the r1 * inverse square factor due to r2
    r = r2 - r1
    return r / (pow(r.Magnitude(), 3) + softeningRadius)


def FastAcceleration(bodies, t, softeningRadius):
    # calculates acceleration on each body. Optimised by not reusing opposite sign pairs of accelerations

    t -= 1
    n = len(bodies)

    # calculate 1/r^2 on each the acceleration of each body for elements that will not be overwritten later
    m = 0
    for i in range(n - 1):
        for j in range(m, n - 1):
            bodies[i].acceleration[j] = InverseSq(bodies[i].position[t], bodies[j + 1].position[t], softeningRadius)

        m += 1

    # reuse the 1/r^2's calculated before for pairs and multiply by secondary mass
    m = 1
    for i in range(1, n):
        for j in range(m):
            bodies[i].acceleration[j] = bodies[j].acceleration[i-1] * -1 * bodies[j].mass

        m += 1

    # multiple 1/r^2 by the secondary mass in the remaining elements
    m = 0
    for i in range(n - 1):
        for j in range(m, n - 1):
            bodies[i].acceleration[j] = bodies[i].acceleration[j] * bodies[j + 1].mass

        m += 1

    #sum for the total acceleration and multiply by G
    for i in range(n):
        bodies[i].totalAcceleration = bodies[i].acceleration.Sum() * G


def Acceleration(bodies, t, softeningRadius):
    # Calculate the acceleration of each body. Not optimised
    t -= 1

    for i in range(len(bodies)):
        #reset total acceleration for the body
        acceleration = Vector3.Zero()

        #add the acceleration component for each secondary body
        for j in range(len(bodies)):
            if i != j:
                acceleration += InverseSq(bodies[i].position[t], bodies[j].position[t], softeningRadius) * bodies[j].mass

        bodies[i].totalAcceleration = acceleration * G


def FirstPrinciples(bodies, t, deltaTime, softeningRadius):
    #Calculate the new position and velocity at time step t
    FastAcceleration(bodies, t, softeningRadius)

    #Set new positions and velocities
    for body in bodies:
        body.velocity[t] = body.velocity[t-1] + body.totalAcceleration * deltaTime
        body.position[t] = body.position[t-1] + body.velocity[t-1] * deltaTime + body.totalAcceleration * pow(deltaTime, 2) * 0.5


def FAcceleration(x, bodies, t, i, softeningRadius):
    # Calculate the acceleration used in F()
    acceleration = Vector3.Zero()

    for j in range(len(bodies)):
        if i != j:
            acceleration += InverseSq(x, bodies[j].position[t-1], softeningRadius) * bodies[j].mass

    return acceleration * G


def F(vectorArray, bodies, t, i, softeningRadius):
    return Vector3Array([vectorArray[1], FAcceleration(vectorArray[0], bodies, t, i, softeningRadius)])


def Euler(bodies, t, deltaTime, softeningRadius):
    FastAcceleration(bodies, t, softeningRadius)

    for body in bodies:
        body.velocity[t] = body.velocity[t - 1] + body.totalAcceleration * deltaTime
        body.position[t] = body.position[t - 1] + body.velocity[t-1] * deltaTime


def SemiEuler(bodies, t, deltaTime, softeningRadius):
    # A more accurate version of the Euler method


    FastAcceleration(bodies, t, softeningRadius)

    for body in bodies:
        body.velocity[t] = body.velocity[t - 1] + body.totalAcceleration * deltaTime
        body.position[t] = body.position[t - 1] + body.velocity[t] * deltaTime


def RK4(bodies, t, deltaTime, softeningRadius):
    #4th order Ruge-Kutta method


    for i in range(len(bodies)):
        x = Vector3Array([bodies[i].position[t-1], bodies[i].velocity[t-1]])

        k1 = F(x, bodies, t, i, softeningRadius) * deltaTime
        k2 = F(x + k1 / 2, bodies, t, i, softeningRadius) * deltaTime
        k3 = F(x + k2 / 2, bodies, t, i, softeningRadius) * deltaTime
        k4 = F(x + k3, bodies, t, i, softeningRadius) * deltaTime

        bodies[i].position[t] = x[0] + k1[0]/6 + k2[0]/3 + k3[0]/3 + k4[0]/6
        bodies[i].velocity[t] = x[1] + k1[1]/6 + k2[1]/3 + k3[1]/3 + k4[1]/6



def ZeroMomentum(bodies):
    # Transform the velocities so that the system has zero momentum
    n = len(bodies)

    totalMomentum = Vector3.Zero()
    for i in range(n):
        totalMomentum += bodies[i].velocity * bodies[i].mass

    velocityCorrection = totalMomentum / sum([bodies[i].mass for i in range(n)])

    for i in range(n):
        bodies[i].velocity -= velocityCorrection


def CentreOfMass(bodies):
    #Get the centre of mass position
    cM = Vector3.Zero()
    massSum = 0
    for body in bodies:
        cM += body.position * body.mass
        massSum += body.mass

    cM /= massSum

    #transform each position to the centre of mass frame
    for body in bodies:
        body.position -= cM


def InitZeroArrays(bodies, numSteps):
    #initialises zero position and velocity arrays for bodies
    for i in range(len(bodies)):
        bodies[i].ZerosPosVel(numSteps)


def InitSmartAcceleration(bodies):
    #initialises the fast acceleration method


    n = len(bodies) - 1

    zeros = Vector3Array([0 for _ in range(n)])

    for body in bodies:
        body.acceleration = deepcopy(zeros)





def VerticalArray(bodies, unit=1, getPositions=True, getVelocities=True):
    # return vertical vector array, [[x1, x2, ...], [y1, y2, ...], [z1, z2, ...]]
    # divide new array by factor to change units
    positions = np.zeros([len(bodies), 3, len(bodies[0].position)])
    velocities = np.zeros([len(bodies), 3, len(bodies[0].velocity)])

    if getPositions and getVelocities:
        for i in range(len(bodies)):
            positions[i][0] = [bodies[i].position[j].x for j in range(len(bodies[i].position))]
            positions[i][1] = [bodies[i].position[j].y for j in range(len(bodies[i].position))]
            positions[i][2] = [bodies[i].position[j].z for j in range(len(bodies[i].position))]

        for i in range(len(bodies)):
            velocities[i][0] = [bodies[i].velocity[j].x for j in range(len(bodies[i].velocity))]
            velocities[i][1] = [bodies[i].velocity[j].y for j in range(len(bodies[i].velocity))]
            velocities[i][2] = [bodies[i].velocity[j].z for j in range(len(bodies[i].velocity))]

        positions /= unit
        velocities /= unit

        return positions, velocities

    elif getPositions:
        for i in range(len(bodies)):
            positions[i][0] = [bodies[i].position[j].x for j in range(len(bodies[i].position))]
            positions[i][1] = [bodies[i].position[j].y for j in range(len(bodies[i].position))]
            positions[i][2] = [bodies[i].position[j].z for j in range(len(bodies[i].position))]

        positions /= unit

        return positions

    elif getVelocities:
        for i in range(len(bodies)):
            velocities[i][0] = [bodies[i].velocity[j].x for j in range(len(bodies[i].velocity))]
            velocities[i][1] = [bodies[i].velocity[j].y for j in range(len(bodies[i].velocity))]
            velocities[i][2] = [bodies[i].velocity[j].z for j in range(len(bodies[i].velocity))]

        velocities /= unit

        return velocities


def CalcKE(bodies):
    # returns the kinetic energy of each body
    kE = np.zeros((len(bodies), len(bodies[0].velocity)))

    for i in range(len(bodies)):
        for j in range(len(bodies[i].velocity)):
            kE[i][j] = 0.5 * bodies[i].mass * pow(bodies[i].velocity[j].Magnitude(), 2)

    return kE


def CalcGPE(bodies):
    #returns the graviational potential energy of each body


    gPE = np.zeros((len(bodies), len(bodies[0].position)))

    for i in range(len(bodies)):
        for t in range(len(bodies[i].position)):
            for j in range(len(bodies)):
                if i != j:
                    gPE[i][t] += bodies[j].mass / (bodies[i].position[t] - bodies[j].position[t]).Magnitude()

            gPE[i][t] *= -bodies[i].mass * G

    return gPE


def TotalEnergy(bodies):
    #returns the total energy of each body


    kE = CalcKE(bodies)
    gPE = CalcGPE(bodies)

    totalEnergy = kE + gPE

    return totalEnergy