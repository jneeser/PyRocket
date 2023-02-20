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
from matplotlib import pyplot as plt


class ChamberGeometry():
    """[Summary]
    Class to calculate the combustion chamber geometry based on a Rao nozzle. 
    Only full nozzle geometries are currently implemented (no truncated nozzles). 
    Contour can be polttend and saved (recommended).
    """
    def __init__(self, D_c, D_t, D_e, L_cyl, r_1, r_2, r_n, phi_conv, phi_div, phi_e, step_size=0.001):
        self.D_c = D_c                                              # chamber diameter [m]
        self.D_t = D_t                                              # throat diameter [m]
        self.D_e = D_e                                              # exit diamter [m]
        self.L_cyl = L_cyl                                          # cylindircal chamber length [m]
        self.r_1 = r_1                                              # converging section inlet radius [m]
        self.r_2 = r_2                                              # beginning of throat radius  [m]
        self.r_n = r_n                                              # end of throat radius [m]
        self.phi_conv = phi_conv * np.pi / 180                      # convergence angle [deg]    
        self.phi_div = phi_div * np.pi / 180                        # divergence angle [deg]
        self.phi_e = phi_e * np.pi / 180                            # exit angle of parabolic nozzle
        self.expansion_ratio = D_e**2 / D_t**2                      # expansion ratio (Ae/At)
        self.dx = step_size                                         # step size for the discretisation [-]

    def contour(self):
        # generate chamber countour x, y coordiantes
        
        # Cylindrical section
        x_cyl = np.arange(0, self.L_cyl, self.dx)
        y_cyl = np.ones(len(x_cyl)) * self.D_c / 2

        # Converging section radius (at least 3 points)
        n1 = max(round((self.r_1 * self.phi_conv ) / self.dx), 3)     # number of discretisations in the radius
        t = np.linspace(0, self.phi_conv, n1)

        x_B = self.L_cyl + self.r_1 * np.sin(t)
        y_B = (self.D_c / 2 - self.r_1) + self.r_1 * np.cos(t)

        L_r1 = self.r_1 * np.sin(self.phi_conv)

        # Converging section 
        a_C = np.tan(self.phi_conv)
        L_con = ((self.D_c / 2 - self.r_1 + self.r_1 * np.cos(self.phi_conv)) - (self.D_t / 2 + self.r_2 - self.r_2 * np.cos(-self.phi_conv))) / a_C
        x_C = np.arange((self.L_cyl + L_r1 + self.dx), (self.L_cyl + L_r1 + L_con), self.dx)
        y_C = (self.D_c / 2 - self.r_1 + self.r_1 * np.cos(self.phi_conv)) - a_C * np.arange(self.dx, L_con, self.dx);


        # Converging radius up to throat (at least 3 points)
        n2 = max(round((self.r_2 * self.phi_conv) / self.dx), 3)      # number of discretisations in the radius
        t = np.linspace(- self.phi_conv, 0, n2)

        L_r21 = self.r_2 * np.sin(self.phi_conv)
        x_D = self.L_cyl + L_r1 + L_con + L_r21 + self.r_2 * np.sin(t)
        y_D = self.D_t / 2 + self.r_2 - self.r_2 * np.cos(t)


        # Throat (at least 3 points)
        n3  = max(round((self.r_n * self.phi_div) / self.dx), 3)      # number of discretisations in the radius
        t = np.linspace(self.phi_div / n3, self.phi_div, n3)

        L_r22 = self.r_n * np.sin(self.phi_div)

        x_E = self.L_cyl + L_r1 + L_con + L_r21 + self.r_2 * np.sin(t)
        y_E = self.D_t / 2 + self.r_n - self.r_n * np.cos(t)

        # parabolic Rao nozzle 
        # starting points of this section
        P1x = x_E[-1]
        P1y = y_E[-1]
        # parabola factors 
        
        a2  = (np.tan(self.phi_e)**(-1) - np.tan(self.phi_div)**(-1)*self.D_e*0.5 / P1y) * (1 - self.D_e*0.5 / P1y)**(-1)
        a1  = (np.tan(self.phi_div)**(-1) - a2) / (2 * P1y)
        a3  = P1x - a1 * P1y**2 - a2 * P1y
        P2x = a1 * (self.D_e*0.5)**2 + a2 * (self.D_e*0.5) + a3

        # number of points along exit parabola (at least 3)
        n4 = max(round((P2x - P1x)/self.dx), 3)
        
        y_F = np.linspace(P1y, self.D_e/2, n4)
        x_F = a1 * y_F**2 + a2 * y_F + a3
 
        self.x = np.concatenate((x_cyl, x_B, x_C, x_D, x_E, x_F), axis=0)
        self.y = np.concatenate((y_cyl, y_B, y_C, y_D, y_E, y_F), axis=0)

        # create 2 dimensional array containing all coordinates
        xy = [[self.x[i], self.y[i]] for i in range(len(self.x))]
        self.geometry = np.vstack(xy)


    def plot_contour(self, path=False):
        plt.plot(self.x, self.y, color='b')
        plt.plot(self.x, -self.y, color='b')
        plt.xlabel('x_coordinate [m]')
        plt.ylabel('y_coordinate [m]')
        plt.title('Chamber Contour')
        plt.axis('equal')
        if path != False:
            plt.savefig(path + '/chamber_contour.png')
            plt.close()
        else:
            plt.show()



class CoolingGeometry():
    """[Summary]
    Class to calculate the cooling channel shapes along the chamber contour. Currently limited to rectangular channels.
    The input parameters 'h_c', 'psi' and 't_w_i' can be input as functions of the axial chamber coordinate. 
    The location of thermocouples can be set to log the temperature of the 2D thermal sim at certain chamber locations. 
    This is meant for program validation purposes and can be disabled in config.py.
    """
    def __init__(self, chamber_geometry, h_c, psi, t_w_i, t_w_o, n_channels):
        self.geom       = chamber_geometry                                 # points along the chamber contour [x,y]
        self.n_channels = n_channels                                       # number of cooling channels
        
        # some inputs can be functions of x coordinates 
        if type(h_c) == float:
            # check if input is a single float
            self.h_c = np.ones(len(self.geom[:,1])) * h_c               # height of the cooling channels [m]
        elif callable(h_c):
            # check if h_c is a function pointer
            self.h_c = h_c(self.geom[:,0])                              # height of cooling channel as function of x 
        else: 
            raise ValueError('h_c must be a float or a function of x coordinate')
        
        if type(psi) == float:
            # check if phi is a single float
            self.psi = np.ones(len(self.geom[:,1])) * psi                # ratio of coverage of the cooling channels to solid wall
        elif callable(psi):
            # check if psi is a function pointer
            self.psi = psi(self.geom[:,0])                               # psi as a fucntion of chamber length if it is indeed  a fucntion
        else: 
            raise ValueError('psi must be a float or a function of x coordinate')

        if type(t_w_i) == float:
            # check if inner wall thickness is a single float
            self.t_w_i = np.ones(len(self.geom[:,1])) * t_w_i              #  inner chamber wall thickness [m]  
        elif callable(t_w_i):
            # check if psi is a function pointer
            self.t_w_i = t_w_i(self.geom[:,0])                               # t_w_i as a fucntion of chamber length if it is indeed  a fucntion
        else: 
            raise ValueError('t_w_i must be a float or a function of x coordinate')
        
        self.t_w_o = np.ones(len(self.geom[:,1])) * t_w_o              # outer chamber wall thickness [m] 

        self.psi_c = 2 * np.pi * self.psi / n_channels                 # radial section of single cooling channel 
        self.psi_w = 2 * np.pi * (1 - self.psi) / n_channels           # raidal section of single uncooled wall segment


    def channel_geometry(self):
        # generate the cooling channel geometry; area, hydraulic diameter and wall thickness
        self.r_i = self.geom[:,1] + self.t_w_i                                      # inner radius of cooling channels 
        self.r_o = self.geom[:,1] + self.t_w_i + self.h_c                           # outer radius of cooling channels 

        self.A_c = self.psi_c / (self.psi_c + self.psi_w) * np.pi * (self.r_o**2 - self.r_i**2)                                                                 # cooling channel area [m^2]
        self.U_c_top = self.psi_c * self.r_o										# circumference of top cooling channels 
        self.U_c_bottom = self.psi_c * self.r_i										# circumference of bottom cooling channels 
        self.U_c_side = self.h_c * 2												# circumference of side cooling channels 
        self.U_c = (self.U_c_top + self.U_c_bottom + self.U_c_side)*self.n_channels # cooling channel circumference [m]
        self.D_h = 4 * self.A_c / self.U_c											# cooling channel hydraulic diameter [m]
		

    def set_thermocouples(self, x, r):
        # set a number of thermocouples for validation purposes
        # loc refers to which boundary on the 2D mesh the thermocouple rests on 
        # input arrays for lnghtwise (x) and radial coordinates (r) 
        # creates an array of dictionaries for all thermocouples in the cooridnate system of the 2D sections
        thermocouples = []

        def closest_node(node):
        # find closest index to thermocouple location
            dist = (self.geom[:,0] - node)**2
            idx = np.argmin(dist)                    # calculates clostes chamber index to TC location
            return idx

        for i in range(len(r)):
            idx = closest_node(x[i])                # axial location index

            section_x = r[i] * np.sin((self.psi_c[idx] + self.psi_w[idx])/2)
            section_y = r[i] * np.cos((self.psi_c[idx] + self.psi_w[idx])/2)
            
            thermocouples.append([idx, section_x, section_y])

        self.thermocouples = np.vstack(thermocouples)

        

if __name__ == "__main__":
    # test functionality of the Chamber and Cooling Geometry Classes

    D_c = 69.05e-3
    D_t = 24.14e-3
    D_e = 50.7e-3
    L_cyl = 74.77e-3
    r_2 = 10e-3
    r_1 = 72e-3
    r_n = 5e-3
    phi_conv = 30
    phi_div = 19
    phi_e  = 14

    cg = ChamberGeometry(D_c, D_t, D_e, L_cyl, r_1, r_2, r_n, phi_conv, phi_div, phi_e, step_size=0.003)
    cg.contour()
    cg.plot_contour()

    n = 8 
    h_c = 3e-3
    phi = 1/3
    t_w_i = 1e-3
    t_w_o = 1e-3

    c = CoolingGeometry(cg.geometry, h_c, phi, t_w_i, t_w_o, n)
    c.channel_geometry()


