import RegenCooling as rc
import PlottingFunctions as pl
import config

# heattransfer sim
sim = rc.HeatTransfer(config.cea, config.gas, config.chamber.geometry, config.material, config.coolant, config.cooling_geom, config.m_dot, config.m_dot_coolant, config.hot_gas_method, config.cooling_method, config.correction, config.eta_combustion)
sim.run_sim()

# plotting
pl.multi_plot(sim.out.T_wall_i, sim.out.T_wall_o, sim.out.q/1e6, sim.out.T_c, 'T_w,i [K]', 'T_w,o [K]', 'q_dot [MW/m^2]', 'T_c [K]')
