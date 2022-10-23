import numpy as np
import scipy.optimize

import RegenCooling as rc
import config

def fun(c):
    global Niter 
    global target_temperature

    sim = rc.HeatTransfer(config.cea, config.gas, config.chamber.geometry, config.material, config.coolant, config.cooling_geom, config.m_dot, config.m_dot_coolant, config.hot_gas_method, config.cooling_method, c, config.eta_combustion)
    sim.run_sim()
    
    print('Iteration:   ', Niter)
    Niter += 1

    return max(sim.out.T_c) - target_temperature

Niter = 1
target_temperature =  450+273                  # target cooling channel exit temperature [K]
res = scipy.optimize.fsolve(fun, 3)
print('Correction factor:   ', res)