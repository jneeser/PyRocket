import RegenCooling as rc
import PlottingFunctions as pl
import config



# heattransfer sim
method = 'cinjarew'
cooling_method = 'gnielinski'
sim = rc.HeatTransfer(config.cea, config.gas, config.chamber.geometry, config.material, config.coolant, config.cooling_geom, config.m_dot, config.m_dot_ox, method, cooling_method, config.eta_combustion)
sim.run_sim()

# plotting
pl.multi_plot(sim.T_wall_arr, sim.q_arr/1e6, sim.T_c_arr, sim.P_c_arr/1e5, 'T_w,i [K]', 'q_dot [MW/m^2]', 'T_c [K]', 'P_c [bar]')
