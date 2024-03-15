#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

import config
import RegenCooling as rc
import PlottingFunctions as pl

# give terminal message about simulation settings
config.settings.output_msg()
config.output.output_msg()

# heat transfer sim
sim = rc.HeatTransfer(
    cea=config.cea,
    gas=config.gas,
    geometry=config.chamber.geometry,
    material=config.material,
    coolant=config.coolant,
    cooling_geometry=config.cooling_geom,
    m_dot=config.m_dot,
    m_dot_coolant=config.m_dot_coolant,
	T_amb = config.ambient_temp,
    output=config.output,
    settings2D=config.settings,
    model=config.hot_gas_method,
    cool_model=config.cooling_method,
    eta_c_star=config.eta_c_star,
    film = config.film,
)

sim.run()

# plotting
plot = pl.Plotting1D(save_path=config.save_path, save=True, show=False)
plot.temperature_plot()
plot.pressure_plot()
plot.heat_transfer_coeff_plot()
plot.heat_flux_plot()
plot.reynolds_plot()
