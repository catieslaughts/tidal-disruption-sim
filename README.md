README.md for project tidal-disruption-sim

This is the final code for a project titled "Modeling Tidal Disruption Events in a Classical Regime," written for a Physical Computing (cosc 89.18) at Dartmouth College in the Fall term of 2020. The starter code was written by Prof. Bo Zhu with help from his TAs. The driver code(s) for the simulation (found in tidal-disruption-sim/proj/final_toolkit) were written by Rory Schadler and Catherine Slaughter. The model writted in a simplified version of that created in Clerici and Gomboc (2020). Our simulation models a tidal disruption event of a star by a black hole in a strictly classical regime. To first order, the simulations look much like similar professional examples. While not entirely physically rigorous, the code creates beautiful simulations what have potential uses in public outreach and teaching.

To run this code, you will need to have CMAKE installed on your machine. Additionally, Linux-based operating systems (including OSX) require dependencies available in freeglut ('brew install glew freeglut', if you have homebrew). From there, you'll need to run the setup script in ./scripts, there is a .bat script for windows machines and a .sh shell script for linux. The code is run by the command: '.\scripts\run_assignment.bat final_toolkit' for windows or './scripts/run_assignment.sh final_toolkit' for linux/OSX

To change the driver being run, go to tidal-disruption-sim/proj/final_toolkit and edit the main.cpp file. Change the value of driver based on the options in the given switch statement. They are as follows:

1. BlackholeCenteredDriver (driver=1)
   * modeling the tidal disruption event of a star about a black hole, visualized in the rest frame of the black hole
   * values like black hole mass (default=300 solar masses) and star mass can be changed in the associated .h file
   * note that the color of the particles making up the star is based on the relative gravitational force felt. Green=lesser net force, red=greater net force.
2. SolarSystemDriver (driver=2)
   * simple model of the solar system
   * in the driver .h file, the values of sun_size_multiplyer and planet_size_multiplyer can be changed for viewing purposes. If both are set to 1, the sizes of the planets and sun are in scale with the distances
3. StarCenteredDriver (driver=3)
   * the same as BlackHoleCenteredDriver, but visualized in the rest frame of the center of mass of the star.
