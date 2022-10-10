# File with unit tests to validate program performance and inputs

import numpy as np
import thermo 
from matplotlib import pyplot as plt

import config


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


if __name__ == "__main__":
    cooling_fluid_test()
    material_property_test()