#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

import numpy as np
import thermo
import os

from IsentropicRelations import Isentropic
from CEAClass import CEA
from GeomClass import ChamberGeometry, CoolingGeometry
from Output import Output1D, Settings2D

import MaterialLib as matlib
import PropLibrary as proplib


# CEA input
fuel            = 'C3H8'                          # choose existing fuel or oxidiser from rocketcea or create new fuel or oxidiser blend in PropLibrary
ox              = 'N2O'
hot_gas_method  = 'cinjarev'                         # use either 'cinjarev', 'standard-bartz' or 'modified-bartz'

# Operating point of combustion chamber 
eta_combustion  = 0.94
MR              = 8
m_dot           = 8e-3                               # total mass flow [kg/s]
m_dot_f         = m_dot / (MR + 1)
m_dot_ox        = m_dot - m_dot_f
Pc              = 10e5                               # Chamber pressure [Pa]

# Coolant input
cooling_fluid   = ['N2O']                    # needs to be list of str
fluid_mass_frac = [1]    			         # mass fractions of the cooling fluid (needs to add up to 1)
m_dot_coolant   = m_dot_ox                 		 	 # mass flow through the cooling channels
inlet_temp      = 293.15                     		 # inlet temperature [K]
inlet_pressure  = 12.5e5                               # Cooling channel inlet pressure [Pa]

cooling_method  = 'gnielinski'                       # use either 'gnielinski', 'dittus-boelter' or 'dittus-boelter-simple', pay attention to suitable Re number range
ambient_temp    = 288.15							 # ambient temperature for radiation boundary condition [K] (must be type float)!!

# Chamber geometry
D_c       = 24e-3                                    # chamber diamter [m]
D_t       = 6.8e-3                                   # throat diameter [m]
D_e       = 32e-3                                   # exit diamter [m]
L_cyl     = 36e-3                                    # cylindrical chamber length [m]
r_1       = 16e-3                                 # converging section radius [m]
r_2       = 3e-3                                     # throat converging section radius [m]
r_n       = 2e-3                                     # throat diverging radius [m]
phi_conv  = 30                                       # convergence angle [deg]
phi_div   = 28                                       # divergence angle [deg]
phi_e     = 15                                       # exit angle [deg]
step_size = 0.002                                   # step size along the chamber contour [m]
material  = matlib.SS14404                           # use entry from MaterialLibrary. Make sure temperature dependent properties are specified

# cooling channel geometry; h_c, psi, t_w_i can be functions of x 
n     = 8                                           # number of cooling channels [int]
h_c   = 3e-3                                         # radial height of cooling channels [m] CAN BE FLOAT OR FUNCTION
# EXAMPLE of function input for psi: psi = lambda x: 1/3 * (1 - 0.0001 * x) 
psi   = 0.4                                         # fill factor of the cooling channels; fraction of the circumferecne covered by the cooling channels (0 - 1) CAN BE FLOAT OR FUNCTION
t_w_i = 2e-3                                         # inner chamber wall thickness [m] CAN BE FLOAT OR FUNCTION
t_w_o = 2e-3                                         # outer chamber wall thickness [m]
start_idx = -1										 # starting index, use 0 for injector side and -1 for nozzle

# 2D Section Simulation settings
cell_size    = 0.1 * t_w_i                           # cell size in 2D section solver, will heavily impact performance
time_step    = 6e2 * (cell_size)**2 / (material.alpha)       # time step in 2D section solver use 2e2 for IN718 and SS14404 and 2e3 for CuCr1Zr
tolerance    = 1e-2                                  # maximum temperature differnece between time steps
max_iter     = 300									 # maximum number of iterations before termination
save_fig     = True                                  # save figures to folder in Output Class
print_result = True                                  # print intermittant maximum temperatures, useful for DEBUGGING
run_time     = 'steady_state'                        # use 'steady_state' as default. use time in [s] if transient solution is desired

# add thermocouple locations for model validation (effecively temperature logging points in the 2D solution)
log_TC = False                                       # temperature at thermocouple locations is to be logged
TC_x   = [14.1e-3, 38.1e-3, 62.1e-3]                 # x coordinate of thermocouples 
TC_r   = [16e-3, 16e-3, 16e-3]                       # radial position of thermocouples               


# save folder path for output class
save_path = "SectionImages"


######################################
# GAS, COOLANT and GEOMETRY SETUP
######################################

# create folder to store sim files
try:
    os.mkdir(save_path)
except:
    pass
    
# generate chamber geometry
chamber = ChamberGeometry(D_c, D_t, D_e, L_cyl, r_1, r_2, r_n, phi_conv, phi_div, phi_e, step_size)
chamber.contour()
chamber.plot_contour(save_path)
geometry = chamber.geometry

# generate cooling channel geometry
cooling_geom = CoolingGeometry(geometry, h_c, psi, t_w_i, t_w_o, n)
cooling_geom.channel_geometry()
cooling_geom.set_thermocouples(TC_x, TC_r)


# generate settings class for 2D thermal sim
settings = Settings2D(cell_size, time_step, tolerance, max_iter, save_fig, print_result, run_time, start_idx, log_TC, cooling_geom.thermocouples)

# generate output class 
output = Output1D(geometry, save_path)

# calculate hot gas properties from NASA CEA
cea = CEA(fuel, ox, Pc) 
cea.metric_cea_output('throat', MR, chamber.expansion_ratio)

# isentropic properties for the hot gas along the chamber contour
gas = Isentropic(Pc, cea.Tc, cea.gamma, cea.Pr, geometry[:,0], geometry[:,1])
gas.calculate()

# coolant Properties from thermo library
coolant = thermo.Mixture(cooling_fluid, ws=fluid_mass_frac, P=inlet_pressure, T=inlet_temp)     