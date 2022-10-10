from scipy import interpolate
import numpy as np
from matplotlib import pyplot as plt


class Material():
	def __init__(self, name, alpha, k, rho, Cp, v, E, a, eps, sig_u):
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
		
	def youngs_modulus(self, E_arr, T_arr):
		self.E = interpolate.interp1d(T_arr, E_arr, kind='linear', fill_value='extrapolate')
		
	def ultimate_strength(self, sig_u_arr, T_arr):
		self.sig_u = interpolate.interp1d(T_arr, sig_u_arr, kind='linear', fill_value='extrapolate')
		
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
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
Ti6Al4V = Material('Ti-6Al-4V', alpha, k, rho, Cp, v, E, a, eps, sig_u)

# set variable material properties from https://www.researchgate.net/publication/299647114_Developments_in_cutting_tool_technology_in_improving_machinability_of_Ti6Al4V_alloy_A_review
Ti6Al4V.youngs_modulus(np.array([114e9, 114e9, 114e9, 85e9, 85e9, 75e9, 42e9]), np.array([293, 373, 473, 573, 673, 773, 873]))
Ti6Al4V.thermal_conductivity(np.array([6.7, 9, 12, 15, 18]), np.array([283, 477, 700, 922, 1144]))
Ti6Al4V.ultimate_strength(np.array([970e6, 900e6, 750e6, 630e6]), np.array([298, 523, 573, 723]))


# IN718 at room temperature
k = 1						# thermal conductivity [W/m/K]
Cp = 1							# specific heat [J/kg/K] at 870 C
rho = 1							# density [kg/m^3]
v = 1                       	# poiosson ratio
E = 1							# youngs modulus [Pa]
a = 1							# thermal expasion coefficinent [1/K]
sig_u = 1                  		# ultimate strength [Pa]
eps = 0.8						# thermal emissivity
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
IN718 = Material('IN718', alpha, k, rho, Cp, v, E, a, eps, sig_u)

# set vairable material properties for IN718
IN718.thermal_conductivity(np.array([11.1, 12.4, 14.12, 16, 17.73, 19.46, 21.19, 23.06]), np.array([294.26, 366.48, 477.59, 588.7, 699.82, 810.92, 922.05, 1033.15]))

# CuCrZr 
k = 1							# thermal conductivity [W/m/K]
Cp = 1							# specific heat [J/kg/K] at 870 C
rho = 1							# density [kg/m^3]
v = 1                       	# poiosson ratio
E = 1							# youngs modulus [Pa]
a = 1							# thermal expasion coefficinent [1/K]
sig_u = 1                  		# ultimate strength [Pa]
eps = 0.8						# thermal emissivity
alpha = k / rho / Cp			# thermal diffusivity [m^2/s]
CuCrZr = Material('CuCrZr', alpha, k, rho, Cp, v, E, a, eps, sig_u)

# set variable material properties for CuCrZr from 'Thermal Fatigue testing of CuCrZr allow for high temperature tooling applications' by Y. Birol
CuCrZr.thermal_conductivity(np.array([290, 315, 315, 310, 290, 270, 310, 330, 350]), np.array([273, 373, 473, 573, 673, 773, 873, 973, 1073]))



