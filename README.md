# PyRocket
PyRocket is a thermal simulaton program for regeneratively cooled rocket engine thrust chambers. It is designed as a reasonably accurate predictive 
tool to assist in design and allow for integration into basic shape optimisation. The program relies on the use of Nusselt number correlations for 
deterimining the heat transfer coefficients from the hot combustion gases, as well as into the cooling fluid. The temperature inside the chamber 
wall is then determined using a 2D transient thermal simulation based on a Finite Volume method. This is repeated for a set number of sections 
along the chamber contour. The predicted heat flux to the cooling fluid in each section is used to calculate a temperature increase of the cooling fluid. 
Basic wall friction estimates are also used to track coolant pressure drop in the channels. 


## Installation
PyRocket relies on several external libraries, as well as GMSH. The installation of dependencies can be done using the install_dependencies.sh file. 


### Installation on WINDOWS 10
The rocektcea library relies on a FORTRAN compiler. It can be installed natively on Windows, but in my 
expereince it is a lot easier to get it to run in a Linux subsystem. 

```bash
    1. as admin, open PowerShell
    2. run: 'wsl --install'
    3. run: 'wsl --install -d Ubuntu-20.04'
    4. run: 'python3 --version' to check the python install (should be above V3.8)

    5. clone PyRocket Files into a folder (locally, not on any network)
    6. open wsl and naviagate to the folder with the PyRocket files 
    7. run: 'sudo bash install_dependencies.sh'
```


### Installation on Ubuntu
The install_dependencies.sh file is set up for installation on Ubuntu, but can easily be altered for use in other Linux distributions. 
```bash
    1. clone PyRocket Files into a folder
    2. open a terminal window and naviagate to the folder with the PyRocket files  
    3. run: 'sudo bash install_dependencies.sh'
```


## Verify Installation
In order to verify that all pacakges were installed correctly, run 'SectionThermalSim.py' in your terminal window. 
```bash
    python3 SectionThermalSim.py
```
This will run a single transient thermal simulation. Do not alter the contents of the config.py file before running. The expected ouput is the maximum temperature
at every time step printed in the terminal window. Additionally a directory called 'TestSectionThermalSim' is created, where an image of the temperature profile 
of the test section is saved. This unit test can be repeated after altering the config.py file for each use case to verify the stability of your new time step and 
mesh resolution settings. 


## Using PyRocket
PyRocket is designed to be used with nearly arbitrary combinations of propellants, cooling fluids and chamber materials. All inputs are changed in config.py 
New propellants can be set in PropLib.py, with new matierals added in MaterialLib.py. The selected coolant must be a fluid present in 
[thermo](https://thermo.readthedocs.io/). A debug mode can be set in config.py to examine the results of the 2D transient simulation to verfiy the chosen mesh and time 
step settings. Run the full simulation by executing 'run.py' in your terminal window. 

```bash
    python3 run.py
```

### Debug mode
Set the 'print_result' option to True in config.py in order to display more detialed information on the each transient thermal simulation. The output will show each time 
step and the resulting maximum temperature. This can show non convergence, in the form of reasonable, but oscillating values, or obvious numerical errors in the 
simulation. In the case of non convergence, the time step can be lowered, or the mesh resolution can be increased. Additionally the maximum iteration count
can be increased if necessary. 


### Setting up chamber geometries
The combustion chamber is defined using a series of parameters in config.py that describe the cylindrical chamber section, and the converging and diverging sections 
respectively. The variable convention here alligns closely with [RPA](https://www.rocket-propulsion.com/index.htm), with the notable differnece, that R1 and R2 are 
exchanged in this program. A plot of the chamber contour is saved in the target directory before each run. The variable 'step_size' in config.py determines the 
number of sections by imposing a minimum resolution in axial direction.  


### Setting up cooling channel geometries
The basic program is set up for a rectangular cooling channel geometry with cooling channel areas as a fraction of local circumference. The basic properties of the cooling 
geometry, such as area and hydraulic diameter are computed in GeomClass.py. The mesh for each chamber section is created in SectionMesh.py, with the final thermal 
simulation and evaluation of results in SectionThermalSim.py. The following naming convention is used for the individual faces (this is relevant for allocating 
boundary conditions and evaluating results).

Boundary Name 		| Description
------------------- | -------------------------------------------------------------
ChamberWall   		| Inner chamber wall in contact with combustion gases
OuterWall     	 	| Outer wall of combustion chamber
CoolantBottomWall	| Wall of the cooling channel closest to the inner chamber wall
CoolantTopWall		| Wall of the cooling channel closest to the outer chamber wall 
CoolantSideWall		| Wall of the cooling channel rib

If the geometry of the cooling channels is fundamentally changed, several functions in several files need to be altered. In order to edit the basic properties, such as 
area and hydraulic diameter, edit the 'CoolingGeometry' class in GeomClass.py. In order to alter the mesh itself, edit the SectionMesh.py file. The 2D mesh geometry is 
produced using [GMSH](https://gmsh.info/). Several tutorials exist on how to use GMSH. Lastly, SectionThermalSim.py needs to be edited. This extends to the faces on which 
boundary conditions are applied, as well as the faces over which boundary fluxes along normal vectors are calculated. These can be found starting in line 52 and 126 respectively. 
Changes to the mesh can be verified by running SectionThermalSim.py


### 2D section solver 
A transient thermal simulation is run on each section of the combustion chamber geomerty. Setting up the finite volume system and boundary conditions is done using 
[FiPy](https://www.ctcms.nist.gov/fipy/). The thermal simualtion takes advantage of symmetry in circumferetial direction, with the flux over the symmetry walls being zero. 

```math
\frac{\partial T}{\partial t} = \alpha \nabla T
```

With 

```math
 q = - k \frac{\partial T}{\partial x}
```

A constant thermal diffusivity is assumed, witht the thermal conductivity dependent on temerpature. The convective heat flux boundary conditions take the form of 

```math
\frac{\partial T}{\partial x} = \frac{1}{k} h_{\alpha} (T_{external} - T_{wall})
```

In varying forms this BC is applied to the boundaries: 'ChmaberWall', 'CoolantBottomWall', 'CoolantTopWall' and 'CoolantSideWall' . The radiation boundary condition
is expressed as

```math
\frac{\partial T}{\partial x} = \frac{\sigma \epsilon}{k} (T_{ambient}^4 - T_{wall}^4)
```

and is applied to the 'OuterWall' boundary. The heat equation is solved for each time step, after which the boundary conditions are updated using the results 
of the previous time step before solving the next time step. The simulation is either terminated when the maximum temperature reaches steady state (within a set tolerance),
or after a set maximum number of iterations is reached. All settings for the 2D thermal sim are found in config.py. During a full run of the program, the boundary flux over
all boundaries is calculated and used to determine the total heat flux, as well as the temperature increase of the cooling fluid as it passes throug each section. 

## Liquid Film Cooling Model
Optional addition of film cooling to supplement regenerative cooling. This is an implementation of the analytical 0D liquid film cooling model proposed by [Shine et al](https://www.researchgate.net/publication/256718300_A_new_generalised_model_for_liquid_film_cooling_in_rocket_combustion_chambers).The model relies on an energy balance between the convective and radiative heat transfer of the combustion gases and the evaporation of the coolant. The model produces an estimate of the liquid film length and assumes that the temperature rise of the wall under the liquid film is very small. The liquid film cooling length and the resultant reduction in convective heat transfer coefficient are passed to the regenerative cooling method. The film cooling is only assumed to be effective as long as the chamber section is within the liquid film length. 

### Plotting Options 
Results of the full chamber simulation are saved in the target directory as a csv file. A select set of saved data is automatically plotted by calling 'PlottingFunctions.py'.
These include the heat transfer coefficients, coolant properties, as well as maximum wall temperature. 


### Chamber Side Heat Transfer Coefficient Models
The heat flux to the combustion chamber is dominated by convective heat transfer, with some contribution from radiation. Three Nusselt number correlations can be chosen for 
the convective heat transfer coefficients. 

* 'cinjarev'
* 'standard-bartz'
* 'modified-bartz'
		
Especially for small chamber geomerties, the 'cinjarev' correlation is recommended.


### Coolant Side Heat Transfer Coefficient Models
Similar to the hot gas side, the convection of heat to the cooling fluid can be estimated using one of three Nusselt number correlations. 

* 'gnielinski'
* 'dittus-boelter'
* 'dittus-boelter-simple'

These correlations have different ranges of applicability in terms of Reynolds number. Check the comments in RegenCooling.py for details. The 'gnielinski' model can lead to
negative Nusselt numbers when the Reynolds number of the cooling fluid is too low. If this error is encountered, use the 'dittus-boelter-simple' correlation instead. The 
'dittus-boelter' correlation also features a viscosity correction for the near-wall fluid film. In practive the 'gnielinski' model is the most accurate, if applicable, as it 
includes information about surface roughness, especially relevant for 3D printed designs. 


### Using Material Library
Four common rocket engine chamber matierals can be selected from MaterialLib.py, 'Ti-6Al-4V', 'IN718', 'CuCr1Zr' and 'SS 1.4404'. These can be selected in config.py.
If a differnet material is desired, it can be added in MaterialLib.py, similar to the existing examples. 


### Creating New Propellant Blends
In config.py any fuel and oxidiser present in [RocketCEA](https://rocketcea.readthedocs.io/en/latest/) can be input. If a new propellant is desired, it can be 
created as a string input, as shown in the PropLibrary.py file. Several examples are already present of new chemicals, as well as propellant blends. 


## Returning cooling channels 
Some regeneratively cooled chamber designs feature cooling channles that first move along the chamber in one direction and return through another set of channels
along the coutour. In order to implement such a design feature in PyRocket, the simulation needs to be run twice. First, set the cooling geomety, with the full
number of channels (irrespecive of their direction) and double the coolant mass flow. This ensures both the appropriate coverage of the cooling channels, as well as 
the correct Reynolds number in the channels. By setting 'start_idx' in config.py to 0 the flow of coolant starts at the injector. Run this simulation. 
Using the results from the first run, 'start_idx' is set back to -1 (to start coolant flow at the nozzle) and the starting temperature and pressure of the coolant 
is set to the final values of the first simulation. The simulation can then be run again for the returning channels. The two simulation runs then provide a reasonably accurate 
estimates of the temperatures around both sets of cooling channels. 


## Optimisation and Parametric Sweeps 
PyRocket can be used to optimise a cooling channel geometry, as well as determining the impact of certain parameters on cooling performance. A sequential optimisation 
routine can be implemented using scipy.optimize.minimize. Using a Monte Carlo approach or parametric sweeps can be run in parallel using python multiprocessing. 
This approach is recommended in most cases as single case run times typically exceed 5 mins and are intrinsically difficult to parallelise. 

## Program Validation
PyRocket has been validated using data gathered on 22 N N2O/Hydrocarbon thrusters. Using data from nine steady state tests, an average total heat flux 
error of 8.3% with a standard deviation of 5.6% was determined. The wall temperature has been validated at discrete locations with a mean error of 11.2% and a standard deviation of 10%.
For Validation purposes the axial and radial coordinate of thermocouples along the camber can be set in config.py. The temperature at that location will then be displayed in the 
terminal output for the simulation. Set 'log_TC' to True to enable the output.  


## Common Issues

### Single section takes more than 10 seconds to solve

This can either be caused by the solver not reaching a steady state solution or the time step being set too small. Abort the simulation and try again with a different adaptive time step configuration. The function 'adaptive_time_step' in 'SectionThermalSim.py' is responsible for setting the time step based on the previous solution. For better convergence behaviour reduce the value of 'step_up' in the 'adaptive_time_step' function (default is 1.2). 

### Fluid propery returned as 'None' type

The library [thermo](https://thermo.readthedocs.io/) does not produce reliable results for Prandtl Number, Cp or viscosity near the supercritical region, 
especially for mixtures. The error will occur in the function 'heat_transfer_coeff_coolant' in the file 'RegenCooling.py'. If this issue occurs, use a 
simplified single componant fluid as coolant. You can also write your own coolant class based on the outputs of the thermo Mixture object and replace the
standard implementation. For this changes need to be made in 'RegenCooling.py' and 'config.py'


    
