Welcome to Raf's Strike point simulator!

This repo contains the code made by me (Raf) to simulate the position of the strike points in a tokamak divertor for the group project of the TU/e fusion course 36 FBA.

The repo contains three files: FieldSimulation.py, ReducedFieldSimulation.py, and StrikePointSim.m.
FieldSimulation.py is the original version of the simulation, primarily used during development and for the creation of a video of a time-varying field.
ReducedFieldSimulation.py is a reduced version of FieldSimulation.py, that creates a separately callable function that determines the magnetic field lines.
ReducedFieldSimulation.py is made to be integrable into a larger system, in contrast to FieldSimulation.py, which does not have this ability.
StrikePointSim.m is a translation of ReducedFieldSimulation.py from Python into MatLab, such that the simulation can be used as the plant of a controller in Simulink.
Whereas FieldSimulation.py and ReducedFieldSimulation.py are completely my own original work, the translation between ReducedFieldSimulation.py and StrikePointSim.m was done using Claude AI.
The validity of the translation was tested by comparing the results of the simulations for equal input parameters. They were found to be equivalent.
