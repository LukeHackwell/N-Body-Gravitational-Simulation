from copy import deepcopy
import numpy as np

G = 6.6743 * pow(10, -11)


def Zeros(n):
    # creates an array of zeros of length n
    return [0 for _ in range(n)]


class Vector3:
    #stores 3D vectors

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, vector2):
        # returns self + vector2
        return Vector3(self.x + vector2.x, self.y + vector2.y, self.z + vector2.z)

    def __sub__(self, vector2):
        # returns subtracts vector2 from self
        return Vector3(self.x - vector2.x, self.y - vector2.y, self.z - vector2.z)

    def __mul__(self, factor):
        # returns self multiplied by scalar factor
        return Vector3(self.x * factor, self.y * factor, self.z * factor)

    def __truediv__(self, factor):
        # returns self divided by scalar factor
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
        # return the magnitude of self
        return pow(pow(self.x, 2) + pow(self.y, 2) + pow(self.z, 2), 0.5)

    @staticmethod
    def Zero():
        # return a zero vector
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


def InverseSq(r1, r2, softeningRadius):
    # returns the r1 * inverse square factor due to r2
    r = r2 - r1
    return r / (pow(r.Magnitude(), 3) + softeningRadius)



def FastAcceleration(bodies, softeningRadius):
    # calculates acceleration on each body. Optimised by not reusing opposite sign pairs of accelerations

    n = len(bodies)

    # calculate 1/r^2 on each the acceleration of each body for elements that will not be overwritten later
    m = 0
    for i in range(n - 1):
        for j in range(m, n - 1):
            bodies[i].acceleration[j] = InverseSq(bodies[i].position, bodies[j + 1].position, softeningRadius)

        m += 1

    # reuse the 1/r^2's calculated before for pairs and multiply by secondary mass
    m = 1
    for i in range(1, n):
        for j in range(m):
            bodies[i].acceleration[j] = bodies[j].acceleration[i - 1] * -1 * bodies[j].mass

        m += 1

    # multiple 1/r^2 by the secondary mass in the remaining elements
    m = 0
    for i in range(n - 1):
        for j in range(m, n - 1):
            bodies[i].acceleration[j] = bodies[i].acceleration[j] * bodies[j + 1].mass

        m += 1

    # sum for the total acceleration and multiply by G
    for i in range(n):
        bodies[i].totalAcceleration = bodies[i].acceleration.Sum() * G


def Acceleration(bodies, softeningRadius):
    # Calculate the acceleration of each body. Not optimised

    for i in range(len(bodies)):
        # reset total acceleration for the body
        acceleration = Vector3.Zero()

        # add the acceleration component for each secondary body
        for j in range(len(bodies)):
            if i != j:
                acceleration += InverseSq(bodies[i].position, bodies[j].position, softeningRadius) * bodies[j].mass

        bodies[i].totalAcceleration = acceleration * G


def FirstPrinciples(bodies, deltaTime, softeningRadius):
    # Calculate the new position and velocity

    FastAcceleration(bodies, softeningRadius)

    # Set new positions and velocities
    for body in bodies:
        oldVelocity = body.velocity
        oldPosition = body.position
        
        body.velocity = oldVelocity + body.totalAcceleration * deltaTime
        body.position = oldPosition + oldVelocity * deltaTime + body.totalAcceleration * pow(deltaTime, 2) * 0.5


def FAcceleration(x, bodies, i, softeningRadius):
    # Calculate the acceleration used in F()

    acceleration = Vector3.Zero()

    for j in range(len(bodies)):
        if i != j:
            acceleration += InverseSq(x, bodies[j].position, softeningRadius) * bodies[j].mass

    return acceleration * G


def F(vectorArray, bodies, i, softeningRadius):
    return Vector3Array([vectorArray[1], FAcceleration(vectorArray[0], bodies, i, softeningRadius)])


def Euler(bodies, deltaTime, softeningRadius):
    FastAcceleration(bodies, softeningRadius)

    for body in bodies:
        oldPosition = body.position
        oldVelocity = body.velocity
        
        body.velocity = oldVelocity + body.totalAcceleration * deltaTime
        body.position = oldPosition + body.oldVelocity * deltaTime


def SemiEuler(bodies, deltaTime, softeningRadius):
    # A more accurate version of the Euler method

    FastAcceleration(bodies, softeningRadius)

    for body in bodies:
        oldPosition = body.position
        oldVelocity = body.velocity
        
        body.velocity = oldVelocity + body.totalAcceleration * deltaTime
        body.position = oldPosition + body.velocity * deltaTime


def RK4(bodies, deltaTime, softeningRadius):
    #4th order Ruge-Kutta method

    for i in range(len(bodies)):
        oldPosition = bodies[i].position
        oldVelocity = bodies[i].velocity
        
        x = Vector3Array([oldPosition, oldVelocity])

        k1 = F(x, bodies, i, softeningRadius) * deltaTime
        k2 = F(x + k1 / 2, bodies, i, softeningRadius) * deltaTime
        k3 = F(x + k2 / 2, bodies, i, softeningRadius) * deltaTime
        k4 = F(x + k3, bodies, i, softeningRadius) * deltaTime

        bodies[i].position = x[0] + k1[0] / 6 + k2[0] / 3 + k3[0] / 3 + k4[0] / 6
        bodies[i].velocity = x[1] + k1[1] / 6 + k2[1] / 3 + k3[1] / 3 + k4[1] / 6




def ZeroMomentum(bodies):
    # Transform the initial velocities so that the system has zero momentum

    n = len(bodies)

    totalMomentum = Vector3.Zero()
    for i in range(n):
        totalMomentum += bodies[i].velocity * bodies[i].mass

    velocityCorrection = totalMomentum / sum([bodies[i].mass for i in range(n)])

    for i in range(n):
        bodies[i].velocity -= velocityCorrection


def CentreOfMass(bodies):
    #Transform initial positions into the centre of mass frame

    cM = Vector3.Zero()
    massSum = 0
    for body in bodies:
        cM += body.position * body.mass
        massSum += body.mass

    cM /= massSum

    #transform each position to the centre of mass frame
    for body in bodies:
        body.position -= cM



def InitSmartAcceleration(bodies):
    #initialises the fast acceleration method

    n = len(bodies) - 1

    zeros = Vector3Array([0 for _ in range(n)])

    for body in bodies:
        body.acceleration = deepcopy(zeros)


def VerticalArray(bodies, unit=1, getPositions=True, getVelocities=True):
    # return vertical vector array: [[x1, x2, ...], [y1, y2, ...], [z1, z2, ...]]
    # divide new array by unit to change units

    positions = np.zeros([3, len(bodies)])
    velocities = np.zeros([3, len(bodies)])

    if getPositions and getVelocities:
        for i in range(len(bodies)):
            positions[0][i] = bodies[i].position.x
            positions[1][i] = bodies[i].position.y
            positions[2][i] = bodies[i].position.z

        for i in range(len(bodies)):
            velocities[0][i] = bodies[i].velocity.x
            velocities[1][i] = bodies[i].velocity.y
            velocities[2][i] = bodies[i].velocity.z

        return positions, velocities

    elif getPositions:
        for i in range(len(bodies)):
            positions[0][i] = bodies[i].position.x
            positions[1][i] = bodies[i].position.y
            positions[2][i] = bodies[i].position.z

        return positions

    elif getVelocities:
        for i in range(len(bodies)):
            velocities[i][0] = bodies.velocity.x
            velocities[i][1] = bodies.velocity.y
            velocities[i][2] = bodies.velocity.z

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
            
            

