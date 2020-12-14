# PBS HS2020 - Project "PBF Splashy Waves"

This is the project for master course Physically-Based Simulation in Computer Graphics offered at ETHZ in 2020 autumn semester. We are group 17 and we implement water splashing/sloshing in container with Position Based Dynamics method (PBD), which resembles Smoothed Particle Hydrodynamic (SPH), but inherits the stability of the geometric, position based dynamics method, allowing large time steps suitable for real-time applications. The paper which we follow can be found here: https://mmacklin.com/pbf_sig_preprint.pdf.

The simulation starts with an amount of fluid particles initialized in the air and then drop into the container. The container moves along x-axis periodically. Both the fluid fall and the movement of container produce fluid splashing/sloshing.

## Installation

### Git and CMAKE
The installation of this project is the same as of other tutorial/homework projects provided in the course. Simply clone the project using Git, and then build the project with CMake with default settings. If you work with Microsoft Visual Studio, please open the PBS.sln file and then set pbf-splashy-waves as start-up project. Then you are free to run the simulation.

### Release and Debug
In order to achieve better visual effect, we recommand to run the simluation in release mode. However, in release mode, the `libigl` UI might not show up correctly. In Debug mode the UI always shows up but the simluation performs significantly slower.

## Parameters
You can modify the following values in FluidSim.h to tune the simulation.

- TIME_STEP_SIZE:               time step between two iterations
- PARTICLE_DISTANCE:            distance between two neighboring fluid particles when initializing
- PARTICLE_RADIUS:              how large each particle should be rendered
- NEIGHBOURHOOD_RADIUS:         smoothing length for SPH kernel function
- BOUNDARY_PARTICLE_DISTANCE:   distance between two neighboring boundary particles when initializing
- PARTICLES_PER_CUBE_SIDE:      number of fluid particles per cube side when intializing
- BOUNDARY_PARTICLE_COLOR:      color of boundary particles
- RENDER_ONLY_FLUID:            whether to only render fluid particles or render both fluid and boundary particles
- ThreadCount:                  degree of parallelization
- MOVING_BOUNDARY:              whether the boundary should move periodically or not move
- AMPLITUDE:                    amplitude of boundary movement
- PERIOD:                       period of boundary movement
- DYNAMIC_PARTICLE_COLORING:    whether the colors of fluid particles are determined by their velocity and y-position at runtime or not 