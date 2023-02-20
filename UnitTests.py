#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

# File with unit tests to validate program performance and inputs

import numpy as np
import thermo 
from matplotlib import pyplot as plt

# functional imports
import config
from Output import Settings2D
from PlottingFunctions import multi_plot

# import functions to be tested
from IsentropicRelations import Isentropic
from SectionThermalSim import HeatEquationSolver



def cooling_fluid_test():
    N = 100                                 # number of steps
    Cp = np.ndarray(N)                      # specific heat at constant pressure [J/kg/K]
    k = np.ndarray(N)                       # thermal conductivity [W/m/K]
    rho = np.ndarray(N)                     # density [kg/m^3] 
    Pr = np.ndarray(N)                      # Prandtl number


    if config.coolant.phase == 'l':
        T_arr = np.linspace(288, 600, N)      # Temperature array [K]
    elif config.coolant.phase == 'g':  
        T_arr = np.linspace(288, 1200, N)  

    for i in range(len(T_arr)):
        # calculate coolant properties at every temperature 
        coolant = thermo.Mixture(IDs=config.coolant.IDs, ws=config.coolant.ws, T=T_arr[i], P=config.coolant.P)
        Cp[i] = coolant.Cp
        rho[i] = coolant.rho
        Pr[i] = coolant.Pr
        k[i] = coolant.Cp * coolant.mu / coolant.Pr

    f, axes = plt.subplots(4, 1)

    axes[0].plot(T_arr, Cp)
    axes[0].set_ylabel('Cp [J/kg/K]')

    axes[1].plot(T_arr, k)
    axes[1].set_ylabel('k [W/m/K]')

    axes[2].plot(T_arr, rho)
    axes[2].set_ylabel('rho [kg/m^3]')

    axes[3].plot(T_arr, Pr)
    axes[3].set_ylabel('Pr')
    
    plt.xlabel('Temperature [K]')
    plt.show()


def material_property_test():
    try:
        config.material.plot(config.material.sig_u, 'ultimate strength [Pa]')
    except:
        print('Warning: no variable ultimate strength of material set')
    try:
        config.material.plot(config.material.E, 'youngs modulus [Pa]')
    except:
        print('Warning: no variable youngs modulus of material set')
    try:
        config.material.plot(config.material.k, 'thermal conductivity [W/m/K]')
    except:
        print('Warning: no variable thremal conductivity of material set')


def isentropic_relations_test():
	# Test isentropic relations class
	gas = Isentropic(config.Pc, config.cea.Tc, config.cea.gamma, config.cea.Pr, config.geometry[:,0], config.geometry[:,1])
	gas.calculate()
	
	multi_plot(gas.M, gas.T_s, gas.p_s/1e5, gas.T_aw, 'M', 'T_s [K]', 'p_s [bar]', 'T_aw [K]')


def section_thermal_sim_test():
    # Test 2D thermal sim 

    def halpha_func(T, idx):
        return 1600, 2400

    def halpha_c_func(T, idx):
        return 1113

    q_rad = 0
    T_c = 288
    T_amb = 288
    idx = 1
    settings = Settings2D(cell_size=2e-4, time_step=4, tolerance=1e-3, max_iter=100, save_fig=False, print_result=False, run_time='steady_state')
    settings.output_msg()
   
	
    solver = HeatEquationSolver(
        idx,
        config.gas,
        config.material,
        config.cooling_geom,
        halpha_func,
        halpha_c_func,
        q_rad,
        T_c,
        T_amb,
        path='TestSectionThermalSim',
        settings=settings
    )
	
    solver.run_sim()
    
    if solver.diff < 1e-3:
        print('Solver converged sucessfully')



if __name__ == "__main__":
    cooling_fluid_test()
    material_property_test()
    isentropic_relations_test()
    section_thermal_sim_test()