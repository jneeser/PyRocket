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
import thermo
import scipy.optimize 

from SectionThermalSim import HeatEquationSolver



class HeatTransfer():
    def __init__(self, cea, gas, geometry, material, coolant, cooling_geometry, m_dot, m_dot_coolant, T_amb, output, settings2D, model='standard-bartz', cool_model='gnielinski', eta_c_star=0.92, film=False):
        """[summary]
        Class to calculate heat tranfer coefficients and radiative heat transfer at arbitrary locations along the chamber contour.
		These functions are passed to the 2D thermal simulation, after which the program advances to the next chamber section. Several options are 
		available for coolant and heat transfer models. Coolant flow properties stored as total conditions, and use a 'thermo' 
		Mixture object to update coolant properties. Hot gas properties from CEA and taken at total chamber conditions. 
		Temperature and pressure of the gas are determined using isentropic expansion.        
        Material thermal condutivity can be temperature dependent.
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
        self.eta_c_star = eta_c_star       			# combustion efficiency
        self.cooling_geometry = cooling_geometry	 	# cooling channel geometry class

        self.A_c = cooling_geometry.A_c					# cooling channel area [m^2]
        self.D_h = cooling_geometry.D_h				    # cooling channel hydraulic diameter [m]

        self.film = film

        self.T_amb = T_amb								# ambient temperature [K], must be float!!!
        self.boltzmann = 5.67e-8						# stefan boltzmann constant

        self.out = output								# output class
        self.settings2D = settings2D					# settings for 2D thermal sim



    def heat_trans_coeff_gas(self, T_wall, idx):
        # get hot gas properties for this chamber location. CEA properties assumed constant, with isentropic gas properties being taken at local point
        gamma = self.cea.gamma
        mu    = self.cea.mu
        cp    = self.cea.Cp
        Pr    = self.cea.Pr
        T_s   = self.gas.T_s[idx]
        p_s   = self.gas.p_s[idx]
        M     = self.gas.M[idx]
        T_aw  = self.gas.T_aw[idx]

        # local geometrical properties 
        local_area  = np.pi * self.geometry[idx,1]**2
        throat_area = np.pi * np.amin(self.geometry[:,1])**2
        D_t         = 2 * np.amin(self.geometry[:,1])
        
        cstar = self.cea.chamber_pressure * throat_area / self.m_dot
        
        # hot gas temperature estimate used for the cinjarev correlation
        T_hg       = T_s + 0.9*(T_aw * self.eta_c_star**2 - T_s)

        # in case of film 
        # stanton number correction without film cooling is 1.0
        St_St0 = 1.0
        if self.film != False:
            if idx in self.film.cooled_idx:
                St_St0 = self.film.St_St0
            else:
                pass
        else:
            pass

        # Nusselt number correlation for the combustion gases
        if self.model == 'standard-bartz':
            # standard bartz correlation, mostly applicable for larger engines. Tends to oveerpredict heat flux near the troat, especially in smaller engines
            sigma  = ((0.5 * T_wall/self.cea.Tc * (1+(gamma-1)/2 * M**2) + 0.5)**(0.68) * (1+(gamma-1)/2 * M**2)**(0.12))**(-1)
            halpha = 0.026/(D_t**0.2) * (mu**0.2)*cp/(Pr**0.6) * (self.gas.p_t/cstar)**0.8 * (throat_area/local_area)**0.9 * sigma

        elif self.model == 'modified-bartz':
            # modified bartz correlation, mostly applicable for larger engines. Tends to oveerpredict heat flux near the troat, especially in smaller engines
            T_f    = 0.5 * T_wall + 0.28 * T_s + 0.22 * T_aw
            G      = self.m_dot/local_area
            halpha = 0.026 * G**0.8/(D_t)**0.2 * mu**0.2*cp/Pr**0.6 * (self.cea.Tc/T_f)**0.68

        elif self.model == 'cinjarev':
            # cinjarev correlation, more useful for medium to small scale motors
            halpha = 0.01975 * self.cea.k**0.18 * (self.m_dot*cp)**0.82 / (2*self.geometry[idx,1])**1.82 * (T_hg/T_wall)**0.35

        else:
            raise ValueError('Invalid heat transfer method. Select: "standard-bartz", "modified-bartz" or "cinjarev"')

        return halpha * St_St0, T_hg
        
        
    def pressure_drop(self, idx):
        # using Darcy Friction Factor to determine pressure drop. Use surface roughness given in Material Class
        surface_roughness = self.material.roughness
        fd = scipy.optimize.fsolve(lambda f: -2 * np.log10(surface_roughness / (3.7 * self.D_h[idx]) + 2.51 / (self.Re * np.sqrt(f))) -1 / np.sqrt(f), 0.001)
        dp = fd * self.section_length[idx] / self.D_h[idx] * 0.5 * self.coolant.rho * self.v_coolant**2 
        
        return dp, fd

        
    def radiation(self, idx):
        # simple radiation model based on h2o and co2 emissions
        R_c = self.geometry[idx,1]					# chamber radius

        p_h2o = self.cea.mole_fractions[1]['H2O'][0] * self.gas.p_s[idx]
        # radiative heat flux of water molecules 
        q_rad = 5.74 * (p_h2o/1e5*R_c)**0.3 * (self.gas.T_s[idx]/100)**3.5

        # check if CO2 is present in the exhaust
        if '*CO2' in self.cea.mole_fractions[1]:
            p_co2 = self.cea.mole_fractions[1]['*CO2'][0] * self.gas.p_s[idx]
            # radiative heat flux of co2 molecules 
            q_rad += 4 * (p_co2/1e5*R_c)**0.3 * (self.gas.T_s[idx]/100)**3.5
    
        return q_rad


    def heat_trans_coeff_coolant(self, T_wall_coolant, idx):
        # coolant flow velocity based on bulk properteis
        self.v_coolant = self.m_dot_coolant / (self.coolant.rho * self.A_c[idx])
        self.Re = self.coolant.rho * self.v_coolant * self.D_h[idx] / self.coolant.mu
		
		# temperature approximation in the near wall fluid film using a logarithmic mean for log temperature profile
        def get_near_wall_fluid():
            T_avg = (T_wall_coolant - self.coolant.T) / np.log(T_wall_coolant / self.coolant.T)
            fluid = thermo.Mixture(IDs=self.coolant.IDs, ws=self.coolant.ws, T=T_avg, P=self.coolant.P)
            return fluid

		# thermodynamic properties of the near wall fluid, Pr implementation of thermo does not work reliably for mixtures near their critical point. Uses gaseouse Pr if liquid Pr returns 'None'
        if self.coolant.phase == 'l':
            Pr = self.coolant.Prl
            Cp = self.coolant.Cpl
            mu = self.coolant.mul

        elif self.coolant.phase == 'g':
            near_wall_coolant = get_near_wall_fluid()
            Pr = near_wall_coolant.Prg
            Cp = near_wall_coolant.Cpg
            mu = near_wall_coolant.mug

        else:
            raise ValueError('Coolant is neither gaseous, nor liquid. Thermo implementation not suitable for supercritical fluids')
        
        # coolant thermal conductivity
        k  = Cp * mu / Pr

        # Nusselt number correlations for the coolant side 
        if self.cool_model == 'dittus-boelter':
            # Dittus-Boelter equation (valid for Re > 1e4 and 0.7 < Pr < 16700)
            near_wall_coolant = get_near_wall_fluid()
            Nu = 0.027 * self.Re**0.8 * Pr**0.33 * (self.coolant.mu / near_wall_coolant.mu)**(0.14)		

        elif self.cool_model == 'dittus-boelter-simple':
            # Dittus-Boelter equation without viscosity correction (valid for Re > 1e4 and 0.7 < Pr < 160)
            Nu = 0.023 * self.Re**0.8 * Pr**0.4

        elif self.cool_model == 'gnielinski':
            # Gnielinski correlation (valid for 3000 < Re < 5e6 and 0.5 < Pr < 2000)
            _, fd = self.pressure_drop(idx)
            Nu = (fd/8) * (self.Re - 1000)*Pr / (1 + 12.7*(fd/8)**0.5 * (Pr**(2/3) - 1))

        else:
            raise ValueError('Invalid heat transfer method. Select: "dittus-boelter", "dittus-boelter-simple" or "gnielinski"')
        
        # convert Nusselt number to heat transfer coefficient
        halpha = Nu * k / self.D_h[idx]

        if halpha < 0:
            raise ValueError('Negative heat transfer coefficient, check applicability of cooling model to Reynolds number range')

        # check if halpha is an array or a float (artifact of the pressure drop dependence or fluid model)
        if isinstance(halpha, float):
            return halpha
        else:
            return halpha[0]


    def section_heat_flux(self, idx):
        self.q_rad = self.radiation(idx)
        
        # solve 2D section temperature profile, passes functions for heat transfer coefficients for temperature dependence of boundary conditions
                
        print('SOLVING Section Number	', len(self.geometry[:,1]) - idx, ' / ', len(self.geometry[:,1]))
        
        solver = HeatEquationSolver(idx, self.gas, self.material, self.cooling_geometry, self.heat_trans_coeff_gas, self.heat_trans_coeff_coolant, self.q_rad, self.coolant.T, self.T_amb, self.out.folder_path, self.settings2D)
        solver.run_sim()

        self.T_wall_i  = solver.inner_wall_temp	    # inner wall temperature taken as maximum temperature of the 2D profile
        self.halpha    = solver.halpha			    # heat transfer coefficients taken as the of the boundary conditions of the 2D solver at the last timestep
        self.halpha_c  = solver.halpha_c_bottom		# average heat transfer coefficient at bottom wall
        self.T_hg      = solver.T_hg

        # total heat flux into the inner chamber wall
        self.q 		   = np.mean(solver.q_chamber)

        # radiation leaving to ambient 
        self.q_rad_out = np.mean(solver.q_outer)
		
        # ensure that FiPy output object is actually a float, multiplied by the number of channels for total area
        Q = float(solver.Q_c) * self.section_length[idx] * self.cooling_geometry.n_channels

 
        # update cooling fluid and pressure. Expected temperature rise and pressure drop 
        dT             = Q / (self.m_dot_coolant * self.coolant.Cp) 
        dp, _          = self.pressure_drop(idx)
        
        # recalcualte cooling fluid properties, such as density, Pr, Cp etc. 
        self.coolant   = thermo.Mixture(IDs=self.coolant.IDs, ws=self.coolant.ws, T=(self.coolant.T+dT), P=(self.coolant.P-dp))
        
        
    def run(self):
        # arrays for section length and area in between the 2D sections solved in the thermal sim 
        self.section_length     = np.ndarray(len(self.geometry[:,1]))

        # determine the section length and inner chamber surface area in each section
        for i in range(len(self.geometry[:,1])):
        # section number 0 at the injector!!!
			
            if self.settings2D.start_idx == -1:
				# for coolant entering at the bottom of the cooling channels, use idx as inverse of i 
                idx = len(self.geometry[:,1]) - i - 1
			
            elif self.settings2D.start_idx == 0:
				# in case coolant starts at injector side
                idx = i
			
            else:
                raise ValueError('Invalid starting point for cooling fluid, select -1 or 0 for nozzle or injector side, respectively')

            if i == 0:
                self.section_length[idx] = 0
                
            else:
                x_len = (self.geometry[idx+1,0] - self.geometry[idx,0])**2
                y_len = (self.geometry[idx+1,1] - self.geometry[idx,1])**2
                
                self.section_length[idx] = np.sqrt(x_len + y_len)

            # solve section heat transfer and temperature field and update coolant properties
            self.section_heat_flux(idx)

            # write to output class
            self.out.halpha[idx]    = self.halpha
            self.out.halpha_c[idx]  = self.halpha_c
            self.out.q_rad[idx]     = self.q_rad
            self.out.q[idx]         = self.q
            self.out.T_wall_i[idx]  = self.T_wall_i
            self.out.T_c[idx]       = self.coolant.T
            self.out.P_c[idx]       = self.coolant.P
            self.out.Re[idx]	    = self.Re
            self.out.T_hg[idx]	    = self.T_hg 
            self.out.v_coolant[idx] = self.v_coolant

        # write to file in output folder
        self.out.write_csv()


