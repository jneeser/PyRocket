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
import scipy.optimize 


class Isentropic():
	"""[Summary]
	Class to calculate isentropic expansion of the combustion gases, as well as the adiabatic wall temperature. 
	((I know T_aw is not an isentropic process)) T_aw includes an assumed recovery factor depending on if the gas
	is in the converging or diverging section. 
	"""
	def __init__(self, p_t, T_t, gamma, Pr, x_coordiantes, y_coordiantes):
		self.p_t = p_t                      # total pressure
		self.T_t = T_t                      # total temperature
		self.gamma = gamma
		self.Pr = Pr
		self.x = x_coordiantes
		self.y = y_coordiantes

		self.M   = np.ndarray(len(self.x))
		self.p_s = np.ndarray(len(self.x))
		self.T_s = np.ndarray(len(self.x))
		self.T_aw = np.ndarray(len(self.x))

	def mach(self):
		A_t = np.amin(self.y)**2 * np.pi			# throat area
		initial_guess = 0.0001
		
		for i in range(len(self.M)):
			A_c = self.y[i]**2 * np.pi

			if A_c == A_t:
				initial_guess = 10

			mach = lambda M:  1/(M*M) * (2/(self.gamma+1) * (1 + (self.gamma-1)/2*M*M))**((self.gamma+1)/(self.gamma-1)) - (A_c/A_t)**2
			self.M[i] = scipy.optimize.fsolve(mach, initial_guess)

	def adiabatic_wall_temp(self):
		# Assumes turbulent flow in the chamber and laminar flow after the throat
		
		for i in range(len(self.T_aw)):
			if self.M[i] >= 1:
				r = self.Pr**(1/2) 
			else:
				r = self.Pr**(1/3) 	
			self.T_aw[i] = self.T_t * (1 + r*(self.gamma - 1)/2 * self.M[i]**2) / (1 + (self.gamma - 1)/2 * self.M[i]**2)

	def pressure(self):
		for i in range(len(self.p_s)):
			self.p_s[i] = self.p_t/((1 + (self.gamma-1)/2 * self.M[i]**2)**(self.gamma/(self.gamma-1)))

	def temperature(self):
		for i in range(len(self.T_s)):
			self.T_s[i] = self.T_t/(1 + (self.gamma-1)/2 * self.M[i]**2)

	def calculate(self):
		self.mach()
		self.pressure()
		self.temperature()
		self.adiabatic_wall_temp()
