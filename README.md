# -MaxwellSolver_NERS570

I propose to revisit the past FDTD program I created by completely recoding it,
and applying the development principles I learned from this class to improve performance, readability,
and alignment with standard development principles. I again will be basing my solutions on the
differential form of Maxwell’s Equations. To avoid working with very large and small numbers, I
elect to write the program primarily using the Gaussian system of units and converting to SI units
only when results are displayed to the user. Additionally, one of my physics professors has agreed to
answer all questions that will arise as I complete this project.

The development process of this project will not start by solving the system in the three-dimensional
case. Instead, the system will first be modeled in two dimensions, after which, the results will be
tested. When, and only when, these tests are passed will development move on to solving the system
in all three dimensions. This progression is designed to streamline the debugging process, allowing for
the refinement of functions and objects in a controlled and manageable environment. For the same
reasons, this program will not immediately support solving physical systems with changing charge
and/or current distributions until results are verified for the static cases.

After this, the implementation of the FDTD algorithm will be standard. A three-dimensional mesh
(A 3D array-like object) will be created which will contain all physical information about the system.
That is to say, each cell of the mesh will contain its own value for the electric and magnetic strength
and direction and information about its source (if one is present), while the index of the cell correlates
to the region of space that it represents. For dynamic simulations, the program will calculate and
apply forces and wave propagation in discrete time steps, which can be adjusted by the user. Finally,
to enhance the realism of the simulation, the boundaries of the computational domain will be equipped
with perfectly-matched layers, or PMLs. These PMLs absorb electromagnetic waves, mimicking the
effect of an open space and preventing reflections that could distort the simulation results.

#################################
TASK LIST:

Workflows: Create a GitHub Repository containing all source code and detailed descriptions
of program functionalities. This will assist in allowing progress to be tracked and program
functionality to be documented. This will also allow for the project to be developed from multiple
workstations.

• Object Oriented Programming and Design: Verbosely comment on each aspect of the program
to ensure high readability and to simplify functionality to be added to it in the future.

• Object Oriented Programming and Design: Create a UML diagram of the project. This will
ensure that the structure of the program is clear from the start.

• Automated Testing Infrastructure: Create unit tests that compare quantities calculated by the
program with those obtained theoretically.

• Automated Testing Infrastructure: Implement a function that calculates the total energy stored
in the simulation to ensure energy is conserved. This is necessary since the program models a
physical system and energy must be conserved throughout the simulation.

• Parallel Computing: Use OpenMP to optimize loops used by the program, allowing for precise
calculations to be run within a shorter time.
