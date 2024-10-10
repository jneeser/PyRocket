#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

from fipy import CellVariable, Gmsh2D, TransientTerm, DiffusionTerm, Viewer
import numpy as np
from matplotlib import pyplot as plt

from Output import Settings2D
from SectionMesh import get_section_mesh


class HeatEquationSolver():
    """[Summary]
    2D transient thermal sim of the chamber wall cross section. Assumes zero flux boundary conditions at the points of symmetry.
    Steady state is reached when the difference between temperatures at subsequent time steps is lower than a tolerance. 
    Settings object is used to pass informatiion about mesh resolution, time step etc. 
    Class is created for each section of the chamber, creates a mesh, sets boundary condtions and solves a 2D heat equation. 
    Information from the previous time step is used for updating the boundary conditions and thermal conductivity of the material.
    The thermal diffusivity is assumed constant with temperature.
    """
    def __init__(self, idx, gas, material, cooling_geometry, halpha_func, halpha_c_func, q_rad, coolant_temperature, ambient_temperature, path, settings):
        self.idx              = idx                      # location along the chamber contour
        self.gas              = gas                      # isentropic gas object
        self.material         = material                 # material object
        self.cooling_geometry = cooling_geometry         # cooling channel geometry object                    
        self.halpha_func      = halpha_func              # hot gas side heat transfer coefficient function pointer
        self.halpha_c_func    = halpha_c_func            # coolant side heat transfer coefficient function pointer
        self.q_rad            = q_rad                    # radiative heat flux [W/m^2]
        self.T_c              = coolant_temperature      # coolant temperature [K]
        self.T_amb            = ambient_temperature      # ambient temperature around the thruster [K]
        self.settings         = settings
        self.path             = path                     # folder path for figures
        self.boltzmann        = 5.67e-8			         # stefan boltzmann constant 


    def boundary_conditions(self, mesh, phi, x, y):
        # determine boundary heat transfer coefficients for coolant and chamber side
        self.halpha, self.T_hg = self.halpha_func(np.mean(phi.faceValue[mesh.physicalFaces["ChamberWall"]]), self.idx)                    # uses maxumum wall temperature to determine halpha
        self.halpha_c_bottom   = self.halpha_c_func(np.mean(phi.faceValue[mesh.physicalFaces["CoolantBottomWall"]]), self.idx)            # uses maximum wall temperature at cooling channel bottom wall for viscosity correction in halpha_c function (if applicable)
        self.halpha_c_top      = self.halpha_c_func(np.mean(phi.faceValue[mesh.physicalFaces["CoolantTopWall"]]), self.idx)            # uses maximum wall temperature at cooling channel bottom wall for viscosity correction in halpha_c function (if applicable)
        self.halpha_c_side     = self.halpha_c_func(np.mean(phi.faceValue[mesh.physicalFaces["CoolantSideWall"]]), self.idx)            # uses maximum wall temperature at cooling channel bottom wall for viscosity correction in halpha_c function (if applicable)
		
        # parameters are multiplied with mesh face normals
        halpha          = self.halpha * mesh.exteriorFaces                   	# heat transfer coefficient of combustion gas
        halpha_c_bottom = self.halpha_c_bottom * mesh.exteriorFaces                   # heat transfer coefficient of the coolant  
        halpha_c_top    = self.halpha_c_top * mesh.exteriorFaces
        halpha_c_side   = self.halpha_c_side * mesh.exteriorFaces
        q_rad           = self.q_rad * mesh.exteriorFaces                      # radiative heat flux
        T_c             = self.T_c * mesh.exteriorFaces                        # coolant temperature
        T_hg            = self.T_hg * mesh.exteriorFaces                       # gas temperature near the inner chamber wall
        T_amb           = self.T_amb * mesh.exteriorFaces                      # ambient temperature around the engine
        k               = self.material.k(phi.faceValue)                       # thermal conductivity of the material (dependent on temperature)

        # Set up heat flux boundary conditions
        # q = -k delta T

        phi.faceGrad.constrain(
            [(-1 / k * halpha_c_top * (phi.faceValue - T_c)) * mesh.faceNormals],
            mesh.physicalFaces["CoolantTopWall"],
        )
        phi.faceGrad.constrain(
            [(-1 / k * halpha_c_side * (phi.faceValue - T_c)) * mesh.faceNormals],
            mesh.physicalFaces["CoolantSideWall"],
        )
        phi.faceGrad.constrain(
            [(-1 / k * halpha_c_bottom * (phi.faceValue - T_c)) * mesh.faceNormals],
            mesh.physicalFaces["CoolantBottomWall"],
        )

        # hot gas side heat transfer on the inner chamber wall, including contribution from radiation
        phi.faceGrad.constrain(
            [(1 / k * (halpha * (T_hg - phi.faceValue) + q_rad)) * mesh.faceNormals],
            mesh.physicalFaces["ChamberWall"],
        )

        # Radiation boundary condition on the outside wall
        phi.faceGrad.constrain(
            [(-self.boltzmann * self.material.eps / k * (phi.faceValue ** 4 - T_amb ** 4)) * mesh.faceNormals],
            mesh.physicalFaces["OuterWall"],
        )


    def adaptive_time_step(self, T_max, min_val=1, low_val=50, max_val=400):
        # adaptive time step based on temperature development of the last two time steps
        delta_T = T_max[-1] - T_max[-2]
        step_down = 2

        if delta_T < 0:
            self.time_step /= step_down
    
        elif delta_T < min_val:
            self.time_step *= 1.0

        elif delta_T < low_val:
            self.time_step *= self.settings.adaptive_up

        elif delta_T > max_val:
            self.time_step /= step_down


    def run_sim(self):
        # generate mesh
        mesh = get_section_mesh(self.settings.cell_size, self.idx)
        # set solution variable and intial guess
        phi = CellVariable(name="Temperature [K]", mesh=mesh, value=float(self.T_amb))

        # create mesh nodes and boundary conditrions
        x, y = mesh.faceCenters
        self.boundary_conditions(mesh, phi, x, y)

        # Equation to solve
        eq = TransientTerm() == DiffusionTerm(coeff=self.material.alpha)

        running = True
        T_max = [self.T_amb, self.T_amb]
        time = 0
        step = 0

        # initalise first time step with value from settings 
        self.time_step = self.settings.time_step

        while running:
            # update boundary conditions
            self.boundary_conditions(mesh, phi, x, y)

            # update time step
            self.adaptive_time_step(T_max)

            # solve equation and update time step
            eq.solve(var=phi, dt=self.time_step)
            time += self.time_step

            # termination critera when steady state is approximately reached
            diff = abs(T_max[-1] - max(phi))
            if diff < self.settings.tolerance:
                running = False 
				
			# terminate at specific run time if it is set
            if self.settings.run_time != 'steady_state':
                if time >= self.settings.run_time:
                    running = False
			
			# terminate if maximum number of iterations is reached 
            if step >= self.settings.max_iter:
                print('No convergence reached after ', self.settings.max_iter, ' iterations. Error of last time step is ', diff)
                print('Increase mesh resolution or decrease time step for better convergence behaviour')
                running = False
				
            if self.settings.print_result:
                print('time: ', round(time,3), ' [s]	T max: ', round(max(phi), 3), ' [K]')

            T_max.append(max(phi))
            step += 1 
            self.diff = diff

        #self.phi = phi
        self.outer_wall_temp = max(phi.faceValue[mesh.physicalFaces["OuterWall"]])
        self.inner_wall_temp = max(phi.faceValue[mesh.physicalFaces["ChamberWall"]])
		
		# boundary fluxes along normal vectors
        dT_dn_c_top    = np.sqrt(phi.faceGrad[0][mesh.physicalFaces["CoolantTopWall"]]**2 + phi.faceGrad[1][mesh.physicalFaces["CoolantTopWall"]]**2)
        dT_dn_c_bottom = np.sqrt(phi.faceGrad[0][mesh.physicalFaces["CoolantBottomWall"]]**2 + phi.faceGrad[1][mesh.physicalFaces["CoolantBottomWall"]]**2)
        dT_dn_c_side   = np.sqrt(phi.faceGrad[0][mesh.physicalFaces["CoolantSideWall"]]**2 + phi.faceGrad[1][mesh.physicalFaces["CoolantSideWall"]]**2)
        dT_dn_chamber  = np.sqrt(phi.faceGrad[0][mesh.physicalFaces["ChamberWall"]]**2 + phi.faceGrad[1][mesh.physicalFaces["ChamberWall"]]**2)
        dT_dn_outer    = np.sqrt(phi.faceGrad[0][mesh.physicalFaces["OuterWall"]]**2 + phi.faceGrad[1][mesh.physicalFaces["OuterWall"]]**2)
		
		# q = k * dT_dn
        self.q_c_top     = self.material.k(phi.faceValue[mesh.physicalFaces["CoolantTopWall"]]) * dT_dn_c_top
        self.q_c_bottom  = self.material.k(phi.faceValue[mesh.physicalFaces["CoolantBottomWall"]]) * dT_dn_c_bottom
        self.q_c_side    = self.material.k(phi.faceValue[mesh.physicalFaces["CoolantSideWall"]]) * dT_dn_c_side
        self.q_chamber   = self.material.k(phi.faceValue[mesh.physicalFaces["ChamberWall"]]) * dT_dn_chamber
        self.q_outer	 = self.material.k(phi.faceValue[mesh.physicalFaces["OuterWall"]]) * dT_dn_outer
		
		# total coolant side heat transferred [W/m]
        self.Q_c = (2 * np.sum([self.q_c_top[i] * self.settings.cell_size for i in range(len(self.q_c_top))])
		          + 2 * np.sum([self.q_c_bottom[i] * self.settings.cell_size for i in range(len(self.q_c_bottom))])
				  + 2 * np.sum([self.q_c_side[i] * self.settings.cell_size for i in range(len(self.q_c_side))]))

        if self.settings.save_fig:
            # save section temperature profile
            viewer = Viewer(vars=phi, datamin=min(phi), datamax=max(phi))
            viewer.plot(filename=(self.path + "/" + str(self.idx)))
            plt.close()

        if self.settings.log_thermocouples:
            # log temperature value at thermocouple locations for validation purposes
            if self.idx in self.settings.thermocouples[:,0]:

                # find which thermocouple is nearest
                TC = np.where(self.idx == self.settings.thermocouples[:,0])[0][0]

                def closest_node(node, nodes):
                    #determine link between TC position and mesh nodes
                    dist = (nodes - node)**2
                    return np.argmin(dist)
                
                # get x position in this section mesh
                x_mesh = closest_node(self.settings.thermocouples[TC, 1], x[mesh.physicalFaces["SideWall"]])
                print('Thermocouple location:    ', self.idx)
                self.T_TC = phi.faceValue[mesh.physicalFaces["SideWall"]][x_mesh]
                print('Thermocouple temperature: ', self.T_TC)
		
		
		
if __name__ == "__main__":

    import config
    import os
	
    try:
        os.mkdir('TestSectionThermalSim')
    except:
        pass
		
    def halpha_func(T, idx):
        return 1600, 2400

    def halpha_c_func(T, idx):
        return 1113

    q_rad = 0
    T_c = 288
    T_amb = 288
    idx = 1
	
    solver = HeatEquationSolver(
        idx,
        config.gas,
        config.material,
        config.cooling_geom,
        halpha_func,
        halpha_c_func,
        q_rad,
        T_c,
        T_amb,
        path='TestSectionThermalSim',
        settings=config.settings
    )

    solver.run_sim()

