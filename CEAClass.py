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

class CEA():
	def __init__(self, fuel, oxidiser, chamber_pressure):
		"""[summary]
		Calcualtes hot gas properties using CEA and converts to SI units. If errors occur with FORTRAN, resart and try again!
		:param chamber_pressure in [Pa]
		:type fuel = string
		:type oxidiser = string
		"""
		self.chamber_pressure = chamber_pressure
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
		self.mole_fractions = self.ispObj.get_SpeciesMoleFractions(Pc=self.imperial_pressure, MR=mixture_ratio, eps=expansion_ratio, frozen=0, frozenAtThroat=0, min_fraction=5e-05)
		self.Cp             = self.Cp * 4186.8																											# coversion to J/kg/K
		self.R              = 8314.46 / self.MW																											# J/kg/K
		self.mu             = self.visc * 0.0001																										# coversion to Pa*s
		self.k              = self.cond * 418.4e-3																										# coversion to W/m/K
		self.Tc             = self.ispObj.get_Tcomb(Pc=self.imperial_pressure, MR=mixture_ratio)*0.555556										        # coversion to K		
