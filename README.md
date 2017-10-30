# Direct Gravitational N-body Simulation #
Repository for developing and parallelizing a direct gravitational N-body simulation during the practical course "Parallele Programmierung" at the University Hamburg, department of Scientific Computing (https://wr.informatik.uni-hamburg.de/teaching/sommersemester_2017/parallele_programmierung). 

Copyright by Nicholas Hickson-Brown and Michael Eidus unless otherwise stated, please refer to the license for this project for more information or the license header of each individual file. Implementation of the Mersenne Twister is provided by Makoto Matsumoto and Takuji Nishimura, please see their implementation for copyright notice.

## Compiling the source code ##
To compile the source code for the computation make sure that the files contained in the __src__ folder are all in the same place and then run the following command: `gcc -o nbody start.c plummer.c mersenne.c hermite.c output.c ediag.c -lm`.

Alternatively you can use the provided __makefile__.

The source code for the visualization does not need to be individually compiled as the folder __Visualization__ contains a zip-File containing the __N Body Visualization 2.0.exe__ for ease of use.

## Running the simulation ##
To run the simulation you need to execute the file __nbody__ after compiling the source code and specify the desired parameters. Available options are:
* _seed_ - specifies the desired seed for the Mersenne Twister (optional)
* _amount_ - defines the desired amount of particles to be generated (required)
* _timestep_ - specifies the timestep to be used for each iteration (required)
* _endtime_ - defines the time at which the simulation terminates (required)

The order of the parameters needs to be: `./nbody [<seed>] <amount> <timestep> <endtime>`.
If no __seed__ is specified, the seed used to initialize the Mersenne Twister is equal to the Unix-Clock at that point.

## Ouput of the simulation ##
During the execution of the simulation a new folder __"run_YEAR_MONTH_DAY_HOURS:MINUTES:SECONDS"__ will be created, which holds all the data produced by the simulation. Files generated are:
* _"log_YEAR_MONTH_DAY_HOURS:MINUTES:SECONDS.txt"_ - contains all important informations about the current run
* _"initial_conditions.csv"_ - contains mass, positions and velocities for all particles at the start of the simulation
* _"energy_diagnostics.csv"_ - contains kinetic, potential and total energy for all iterations, in that order
* _"iteration_X.csv"_ - subsequent iterations which contain the respective mass, positions and velocities for all particles

The order of the particle information within initial_conditions.csv and the iteration_X.csv files is as follows: 

__x-Axis-Position, y-Axis-Position, z-Axis-Position, Mass, x-Axis-Velocity, y-Axis-Velocity, z-Axis-Velocity__

## Visualizing the generated output ##
To visualize the generated data from the simulation make sure that the executable __N Body Visualization 2.0.exe__, the dll's __freetype6.dll__ and __zlib1.dll__, the folder __shaders__, __fonts__ and the folder containing the generated data are all in the same place. 

Rename the folder __"run_YEAR_MONTH_DAY_HOURS:MINUTES:SECONDS"__ to just __"run"__ and start the executable. Two windows should appear, one containing information and another containing the visualization.

By pressing __alt + enter__ the fullscreen mode can be activated and deactivated.

The visualization can be started and stopped by pressing __space__ and reset by pressing __R__.

By pressing the __arrow keys__ zoom and rotation can be controlled.

By pressing __M__ mouselook can be activated and deactivated.

The point size of the particles can be modified by pressing __1__, __2__, __3__.

By pressing __D__ data can be displayed and hidden.

By pressing __F__ circles for distance measuring can be displayed and hidden.

By pressing __NUM_1__, __NUM_4__, __NUM_2__ and __NUM_5__ the aforementioned circles can be changed in size.

By pressing __S__ slowmotion can be activated and deactivated.

By pressing __C__ colored mode can be activated and deactivated.

By pressing __P__ the experimental painting mode can be activated and deactivated.

Press __escape__ to exit the visualization.
