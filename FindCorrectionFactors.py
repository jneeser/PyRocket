import numpy as np
import scipy.optimize

import RegenCooling as rc
import config

def fun(c):
    
    sim = rc.HeatTransfer(config.cea, config.gas, config.chamber.geometry, config.material, config.coolant, config.cooling_geom, config.m_dot, config.m_dot_coolant, config.hot_gas_method, config.cooling_method, c, config.eta_combustion)
    sim.run_sim()
    print('hi')

    return max(sim.out.T_c) - (510+273)

res = scipy.optimize.fsolve(fun, 3)
print(res)