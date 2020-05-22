# N-Body-Gravitational-Simulation
Provides code to run and analyse N-body gravitational simulations. Includes the following features:
* Set a uniform initial positions distribution.
* Set a uniform spherical velocity distribution with magnitudes of the infered virial velocity distribution.
* Set a power law (Salpeter-like) IMF.
* Run multiple simulations for linearly or logarithmically varying time step length, keeping the total simulated time the same. This is useful for analysing how the time step length affects the accuracy of a simulation.
* Code is parallelised using the multiprocessing module.

NBody.py is used to set the initial conditions of the system, the time step length and the total time to be simulated, and to run the simulations. 

LightSim.py contains classes and functions related to the physics of the simulation.

Plotter.py analyses results from simulations.


