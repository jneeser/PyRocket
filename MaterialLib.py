#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

from scipy import interpolate
import numpy as np
from matplotlib import pyplot as plt


class Material():
	def __init__(self, name, alpha, k, rho, Cp, v, E, a, eps, sig_u, roughness):
		self.name =	name			# str of material name
		self.alpha = alpha			# thermal diffusivity [m^2/s]
		self.k = k					# thermal conductivity [W/m/K]
		self.rho = rho				# density [kg/m^3]
		self.Cp = Cp				# specific heat [J/kg/K]
		self.v = v 					# poiosson ratio
		self.E = E					# youngs modulus [Pa]
		self.a = a					# thermal expasion coefficinent [1/K]
		self.eps = eps				# Thermal emissivity 
		self.sig_u = sig_u			# Ultimate strength [Pa]
		self.roughness = roughness	# effective roughness height [m]
		
	def thermal_conductivity(self, k_arr, T_arr):
		self.k = interpolate.interp1d(T_arr, k_arr, kind='linear', fill_value='extrapolate')
		
	def plot(self, material_property, label):
		T_range = np.linspace(288, 1000, 100)
		plt.plot(T_range, material_property(T_range))
		plt.xlabel('temperature [K]')
		plt.ylabel(label)
		plt.title(self.name)
		plt.show()


# Ti-6Al-4V at room temperature
k = 6.7							# thermal conductivity [W/m/K]
Cp = 930						# specific heat [J/kg/K] at 870 C
rho = 4430						# density [kg/m^3]
v = 0.342                       # poiosson ratio
E = 114e9						# youngs modulus [Pa]
a = 9.7e-6						# thermal expasion coefficinent [1/K]
sig_u = 550e6                   # ultimate strength [Pa]
eps = 0.8						# thermal emissivity
roughness = 2e-5				# effective roughness height [m] from Materialise
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
Ti6Al4V = Material('Ti-6Al-4V', alpha, k, rho, Cp, v, E, a, eps, sig_u, roughness)

# set variable material properties from https://www.researchgate.net/publication/299647114_Developments_in_cutting_tool_technology_in_improving_machinability_of_Ti6Al4V_alloy_A_review
Ti6Al4V.thermal_conductivity(np.array([6.7, 9, 12, 15, 18]), np.array([283, 477, 700, 922, 1144]))


# IN718 at room temperature
k = 16						    # thermal conductivity [W/m/K]
Cp = 501					    # specific heat [J/kg/K] at 400 C
rho = 8170				        # density [kg/m^3]
v = 1                       	# poiosson ratio
E = 1							# youngs modulus [Pa]
a = 1							# thermal expasion coefficinent [1/K]
sig_u = 1                  		# ultimate strength [Pa]
eps = 0.8						# thermal emissivity
roughness = 0.7e-5				# effective roughness height [m] from Materialize
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
IN718 = Material('IN718', alpha, k, rho, Cp, v, E, a, eps, sig_u, roughness)

# set vairable material properties for IN718
IN718.thermal_conductivity(np.array([11.9, 13.7, 16.9, 21.7, 25.6, 22.9, 19.1, 17.7]), np.array([22, 233, 448, 657, 866, 1079, 1289, 1500])+273.15)

# CuCr1Zr 
k = 310							# thermal conductivity [W/m/K] at 400 C
Cp = 520					    # specific heat [J/kg/K] at 400 C
rho = 8900						# density [kg/m^3]
v = 1                       	# poiosson ratio
E = 1							# youngs modulus [Pa]
a = 1							# thermal expasion coefficinent [1/K]
sig_u = 1                  		# ultimate strength [Pa]
eps = 0.8						# thermal emissivity
roughness = 2.3e-5				# effective roughness height [m] from SLM Solutions
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
CuCr1Zr = Material('CuCr1Zr', alpha, k, rho, Cp, v, E, a, eps, sig_u, roughness)

# set variable material properties for CuCr1Zr from 'Thermal Fatigue testing of CuCrZr allow for high temperature tooling applications' by Y. Birol
CuCr1Zr.thermal_conductivity(np.array([290, 292, 300, 302, 305, 310, 315, 320, 330]), np.array([273, 373, 473, 573, 673, 773, 873, 973, 1073]))


# SS 1.4404 
k = 23					        # thermal conductivity [W/m/K] at 500 C
Cp = 590					    # specific heat [J/kg/K] at 500 C
rho = 7750						# density [kg/m^3] at 500 C
v = 1                       	# poiosson ratio
E = 1							# youngs modulus [Pa]
a = 1							# thermal expasion coefficinent [1/K]
sig_u = 1                  		# ultimate strength [Pa]
eps = 0.8						# thermal emissivity
roughness = 0.9e-5				# effective roughness height [m] from Materialize
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
SS14404 = Material('SS14404', alpha, k, rho, Cp, v, E, a, eps, sig_u, roughness)
# variable material properties from 'Transient thermo-mechanical modeling of stress evolution and re-melt volume fraction in electron beam additive manufacturing process' by R.K. Adhitan
SS14404.thermal_conductivity(np.array([13.5, 16, 17.5, 19.5, 22, 23.5, 24.5, 25.5, 27, 28, 29]), np.array([273, 373, 473, 573, 673, 773, 873, 973, 1073, 1173, 1273]))


