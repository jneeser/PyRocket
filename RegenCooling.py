import numpy as np
import thermo
import scipy.optimize 

from Output import Output1D




class HeatTransfer():
	def __init__(self, cea, gas, geometry, material, coolant, cooling_geometry, m_dot, m_dot_coolant, model='standard-bartz', cool_model='gnielinski', correction=1, eta_combustion=0.92):
		"""[summary]
		Coolant flow properties stored as total conditions
		Hot gas properties from CEA
		USE SI UNITS (sorry freedom lovers)
		"""
		self.cea = cea                       			# NASA CEA object 
		self.gas = gas                      			# isentropic gas object 
		self.geometry = geometry             			# engine contour
		self.material = material             			# material object
		self.coolant = coolant               			# themro Mixture object
		self.m_dot = m_dot                   			# total mass flow
		self.m_dot_coolant = m_dot_coolant   			# coolant mass flow [kg/s]
		self.model = model                   			# hot gas side heat transfer model
		self.cool_model = cool_model		 			# coolant side heat transfer model
		self.eta_comb = eta_combustion       			# combustion efficiency
		self.cooling_geometry = cooling_geometry	 	# cooling channel geometry class
		self.correction = correction					# correction factor 

		self.A_c = cooling_geometry.A_c					# cooling channel area [m^2]
		self.D_h = cooling_geometry.D_h				    # cooling channel hydraulic diameter [m]
		self.t_w_i = cooling_geometry.t_w_i  	      	# inner wall thickness [m]
		self.t_w_o = cooling_geometry.t_w_o  	      	# inner wall thickness [m]
		
		self.T_amb = 288								# ambient temperature [K]
		self.boltzmann = 5.67e-8						# stefan boltzmann constant 

		# section property arrays
		self.section_length = np.ndarray(len(geometry[:,1]))
		self.section_area = np.ndarray(len(geometry[:,1]))

		# Output arrays
		self.out = Output1D()



	def heat_trans_coeff_gas(self, T_wall, idx):
		gamma = self.cea.gamma
		mu    = self.cea.mu
		cp    = self.cea.Cp
		Pr    = self.cea.Pr
		T_s   = self.gas.T_s[idx]
		p_s   = self.gas.p_s[idx]
		M     = self.gas.M[idx]
		T_aw  = self.gas.T_aw[idx]

		local_area  = np.pi * self.geometry[idx,1]**2
		throat_area = np.pi * np.amin(self.geometry[:,1])**2
		D_t         = 2 * np.amin(self.geometry[:,1])
		
		cstar = self.cea.chamber_pressure * throat_area / self.m_dot

		self.T_hg      = T_s + 0.8*(T_aw * self.eta_comb**2 - T_s)

		if self.model == 'standard-bartz':
			sigma  = ((0.5 * T_wall/self.cea.Tc * (1+(gamma-1)/2 * M**2) + 0.5)**(0.68) * (1+(gamma-1)/2 * M**2)**(0.12))**(-1)
			halpha = 0.026/(D_t**0.2) * (mu**0.2)*cp/(Pr**0.6) * (self.gas.p_t/cstar)**0.8 * (throat_area/local_area)**0.9 * sigma

		elif self.model == 'modified-bartz':
			T_f    = 0.5 * T_wall + 0.28 * T_s + 0.22 * T_aw
			G      = self.m_dot/local_area
			halpha = 0.026 * G**0.8/(D_t)**0.2 * mu**0.2*cp/Pr**0.6 * (self.cea.Tc/T_f)**0.68

		elif self.model == 'cinjarew':
			halpha = 0.01975 * self.cea.k**0.18 * (self.m_dot*cp)**0.82 / (2*self.geometry[idx,1])**1.82 * (self.T_hg/T_wall)**0.35

		else:
			raise ValueError('Invalid heat transfer method. Select: "standard-bartz", "modified-bartz" or "cinjarew"')

		return halpha 
		
		
	def pressure_drop(self, idx, surface_roughness=6e-6):
		fd = scipy.optimize.fsolve(lambda f: -2 * np.log10(surface_roughness / (3.7 * self.D_h[idx]) + 2.51 / (self.Re * np.sqrt(f))) -1 / np.sqrt(f), 0.001)
		dp = fd * self.section_length[idx] / self.D_h[idx] * 0.5 * self.coolant.rho * self.v_coolant**2 
		
		return dp, fd

		
	def radiation(self, idx):
		R_c = self.geometry[idx,1]					# chamber radius

		p_co2 = self.cea.mole_fractions[1]['*CO2'][0] * self.gas.p_s[idx]
		p_h2o = self.cea.mole_fractions[1]['H2O'][0] * self.gas.p_s[idx]

		q_r_co2 = 4 * (p_co2/1e5*R_c)**0.3 * (self.gas.T_s[idx]/100)**3.5
		q_r_h2o = 5.74 * (p_h2o/1e5*R_c)**0.3 * (self.gas.T_s[idx]/100)**3.5

		return q_r_h2o + q_r_co2


	def heat_trans_coeff_coolant(self, T_wall_coolant, idx):
		self.v_coolant = self.m_dot_coolant / (self.coolant.rho * self.A_c[idx])
		
		Pr = self.coolant.Pr
		Cp = self.coolant.Cp
		mu = self.coolant.mu
		
		#Pr = 0.75 + 1.63/np.log(1+Pr/0.0015) 			# turbulent Pr correction
		self.Re = self.coolant.rho * self.v_coolant * self.D_h[idx] / mu
		k  = Cp * mu / Pr

		if self.cool_model == 'dittus-boelter':
			# Dittus-Boelter equation (valid for Re > 1e4 and 0.7 < Pr < 16700)
			near_wall_coolant = thermo.Mixture(IDs=self.coolant.IDs, ws=self.coolant.ws, T=T_wall_coolant, P=self.coolant.P)
			#Nu = 0.027 * self.Re**0.8 * Pr**0.33 * (mu / near_wall_coolant.mul)**(-0.14)		
			Nu = 0.023 * self.Re**0.8 * Pr**0.4

		elif self.cool_model == 'hess-kunz':
			# from AE4S01 (TRP) reader (unknown validity)
			near_wall_coolant = thermo.Mixture(IDs=self.coolant.IDs, ws=self.coolant.ws, T=T_wall_coolant, P=self.coolant.P)
			Nu = 0.0208 * self.Re**0.8 * Pr**0.4 * (1 + 0.01457 * near_wall_coolant.mu / mu)
		
		elif self.cool_model == 'gnielinski':
			# Gnielinski correlation (valid for 3000 < Re < 5e6 and 0.5 < Pr < 2000)
			_, fd = self.pressure_drop(idx)
			Nu = (fd/8) * (self.Re - 1000)*Pr / (1 + 12.7*(fd/8)**0.5 * (Pr**(2/3) - 1))

		else:
			raise ValueError('Invalid heat transfer method. Select: "dittus-boelter", "hess-kunz" or "gnielinski"')

		
		halpha = Nu * k / self.D_h[idx]

		return halpha
		

	def section_solver(self, idx, maxIter=100, tol=1e-6):
		T_wall_guess = 600
		T_wall_coolant = 400
		Niter = 0
		diff = 1

		while diff > tol:
			self.halpha   = self.heat_trans_coeff_gas(T_wall_guess, idx)
			self.halpha_c = self.heat_trans_coeff_coolant(T_wall_coolant, idx)
			self.halpha_c = self.correction * self.cooling_geometry.channel_efficiency(self.material.k(T_wall_guess), self.halpha_c, idx)
			self.q_rad    = self.radiation(idx)
			
			self.q        = (self.gas.T_aw[idx] - self.coolant.T + self.q_rad/self.halpha) / (1/self.halpha + self.t_w_i[idx]/self.material.k(T_wall_guess) + 1/self.halpha_c)
			self.T_wall_i = -((self.q - self.q_rad)/self.halpha - self.gas.T_aw[idx])
			
			self.T_wall_o = -self.q * self.t_w_i[idx]/self.material.k(T_wall_guess) + self.T_wall_i

			diff   = np.max([abs(self.T_wall_i - T_wall_guess), abs(self.T_wall_o - T_wall_coolant)])
			
			Niter += 1

			if Niter > maxIter:
				raise ValueError('Maximum number of iterations reached without convergence')

			T_wall_guess   = self.T_wall_i
			T_wall_coolant = self.T_wall_o

		# update coolant properties
		dT           = self.q * self.section_area[idx] / (self.m_dot_coolant * self.coolant.Cp) 
		dp, _        = self.pressure_drop(idx)
		
		self.coolant = thermo.Mixture(IDs=self.coolant.IDs, ws=self.coolant.ws, T=(self.coolant.T+dT), P=(self.coolant.P-dp))
		
	def run_sim(self):

		for i in range(len(self.geometry[:,1])):
			idx = len(self.geometry[:,1]) - i - 1

			if i == 0:
				self.section_length[idx] = 0
				self.section_area[idx] = 0
				
			else:
				x_len = (self.geometry[idx+1,0] - self.geometry[idx,0])**2
				y_len = (self.geometry[idx+1,1] - self.geometry[idx,1])**2
				
				self.section_length[idx] = np.sqrt(x_len + y_len)
				self.section_area[idx] = np.pi * 2 *self.geometry[idx,1] * self.section_length[idx]
			
			self.section_solver(idx)


			self.out.halpha[idx]    = self.halpha
			self.out.halpha_c[idx]  = self.halpha_c
			self.out.q_rad[idx]     = self.q_rad
			self.out.q[idx]         = self.q
			self.out.T_wall_i[idx]  = self.T_wall_i
			self.out.T_wall_o[idx]  = self.T_wall_o
			self.out.T_c[idx]       = self.coolant.T
			self.out.P_c[idx]       = self.coolant.P
			self.out.Re[idx]	    = self.Re
			self.out.T_hg[idx]	    = self.T_hg 
			self.out.v_coolant[idx] = self.v_coolant

