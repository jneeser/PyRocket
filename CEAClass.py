#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

from rocketcea.cea_obj import CEA_Obj

class BipropCEA():
	def __init__(self, fuel, oxidiser, chamber_pressure, ambient_pressure=1.0):
		"""Calcualtes hot gas properties using CEA and converts to SI units. If errors occur with FORTRAN, resart and try again!

		Args:
			fuel (str): fuel already present in RocketCEA or from PropLibrary.py
			oxidiser (str): oxidiser already present in RocketCEA or from PropLibrary.py
			chamber_pressure (float): combustion chamber total pressure in Pa
			ambient_pressure (float, optional): ambient pressure in Pa. Defaults to 1.0.
		"""
		self.chamber_pressure = chamber_pressure
		self.ambient_pressure = ambient_pressure
		self.imperial_amb_pressure = 0.000145038 * self.ambient_pressure 			#conversion to psia
		self.imperial_pressure = 0.000145038 * self.chamber_pressure 			#conversion to psia
		self.ispObj = CEA_Obj( oxName=oxidiser, fuelName=fuel)
		self.ispObj.get_full_cea_output()

	def chamber_gas_properties(self, mixture_ratio, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Chamber_Transport(Pc=self.imperial_pressure, MR=mixture_ratio, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_Chamber_MolWt_gamma(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)

	def throat_gas_properties(self, mixture_ratio, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Throat_Transport(Pc=self.imperial_pressure, MR=mixture_ratio, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_Throat_MolWt_gamma(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)
	
	def exit_gas_properties(self, mixture_ratio, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Exit_Transport(Pc=self.imperial_pressure, MR=mixture_ratio, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_exit_MolWt_gamma(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)


	def metric_cea_output(self, location, mixture_ratio, expansion_ratio):
		if location == 'chamber':
			self.chamber_gas_properties(mixture_ratio, expansion_ratio)
		elif location == 'throat':
			self.throat_gas_properties(mixture_ratio, expansion_ratio)
		elif location == 'exit':
			self.exit_gas_properties(mixture_ratio, expansion_ratio)
		else:
			raise ValueError('Invalid location, use "chamber," "throat" or "exit"')

		self.isp, self.cstar, _ =  self.ispObj.getFrozen_IvacCstrTc(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)
		
		self.cstar          = self.cstar * 0.3048																								        # coversion to m/s
		self.ispAmb         = self.ispObj.estimate_Ambient_Isp(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio, Pamb=self.imperial_amb_pressure)
		self.mole_fractions = self.ispObj.get_SpeciesMoleFractions(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio, frozen=0, frozenAtThroat=0, min_fraction=5e-05)
		self.Cp             = self.Cp * 4186.8																											# coversion to J/kg/K
		self.R              = 8314.46 / self.MW																											# J/kg/K
		self.mu             = self.visc * 0.0001																										# coversion to Pa*s
		self.k              = self.cond * 418.4e-3																										# coversion to W/m/K
		self.Tc             = self.ispObj.get_Tcomb(Pc=self.imperial_pressure, MR=mixture_ratio)*0.555556
		self.Pc             = self.chamber_pressure												        # coversion to K		
		self.rho 			= self.Pc / (self.R * self.Tc)	
		self.h              = self.ispObj.get_Chamber_H(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio)	


class MonopropCEA():
	def __init__(self, propellant, chamber_pressure, ambient_pressure=1.0):
		"""Calcualtes hot gas properties using CEA and converts to SI units. If errors occur with FORTRAN, resart and try again!

		Args:
			propellant (str): propellant already present in RocketCEA or from PropLibrary.py
			chamber_pressure (float): combustion chamber total pressure in Pa
			ambient_pressure (float, optional): ambient pressure in Pa. Defaults to 1.0.
		"""
		self.chamber_pressure = chamber_pressure
		self.ambient_pressure = ambient_pressure
		self.imperial_pressure = 0.000145038 * self.chamber_pressure 			#conversion to psia\
		self.imperial_amb_pressure = 0.000145038 * self.ambient_pressure 			#conversion to psia
		self.ispObj = CEA_Obj(propName=propellant)
		self.ispObj.get_full_cea_output()

	def chamber_gas_properties(self, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Chamber_Transport(Pc=self.imperial_pressure, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_Chamber_MolWt_gamma(Pc=self.imperial_pressure, eps=expansion_ratio)

	def throat_gas_properties(self, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Throat_Transport(Pc=self.imperial_pressure, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_Throat_MolWt_gamma(Pc=self.imperial_pressure, eps=expansion_ratio)
	
	def exit_gas_properties(self, expansion_ratio):
		self.Cp, self.visc, self.cond, self.Pr = self.ispObj.get_Exit_Transport(Pc=self.imperial_pressure, frozen=1)
		self.MW, self.gamma                    = self.ispObj.get_exit_MolWt_gamma(Pc=self.imperial_pressure, eps=expansion_ratio)


	def metric_cea_output(self, location, expansion_ratio):
		if location == 'chamber':
			self.chamber_gas_properties(expansion_ratio)
		elif location == 'throat':
			self.throat_gas_properties(expansion_ratio)
		elif location == 'exit':
			self.exit_gas_properties(expansion_ratio)
		else:
			raise ValueError('Invalid location, use "chamber," "throat" or "exit"')

		self.isp, self.cstar, _ =  self.ispObj.getFrozen_IvacCstrTc(Pc=self.imperial_pressure, eps=expansion_ratio)
		
		self.cstar          = self.cstar * 0.3048	
		self.ispAmb         = self.ispObj.estimate_Ambient_Isp(Pc=self.imperial_pressure, eps=expansion_ratio, Pamb=self.imperial_amb_pressure)	
		self.mole_fractions = self.ispObj.get_SpeciesMoleFractions(Pc=self.imperial_pressure, eps=expansion_ratio, frozen=0, frozenAtThroat=0, min_fraction=5e-05)
		self.Cp             = self.Cp * 4186.8																											# coversion to J/kg/K
		self.R              = 8314.46 / self.MW																										# J/kg/K
		self.mu             = self.visc * 0.0001																										# coversion to Pa*s
		self.k              = self.cond * 418.4e-3																										# coversion to W/m/K
		self.Tc             = self.ispObj.get_Tcomb(Pc=self.imperial_pressure)*0.555556	                                                                # coversion to K	
		self.Pc             = self.chamber_pressure		
		self.rho 			= self.Pc / (self.R * self.Tc)						        
