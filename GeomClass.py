import numpy as np
from matplotlib import pyplot as plt


class ChamberGeometry():
    def __init__(self, D_c, D_t, D_e, L_cyl, r_1, r_2, phi_conv, phi_div, step_size=0.0004):
        self.D_c = D_c                                              # chamber diameter [m]
        self.D_t = D_t                                              # throat diameter [m]
        self.D_e = D_e                                              # exit diamter [m]
        self.L_cyl = L_cyl                                          # cylindircal chamber length [m]
        self.r_1 = r_1                                              # converging section inlet radius [m]
        self.r_2 = r_2                                              # diverging section radius  [m]
        self.phi_conv = phi_conv * np.pi/180                        # convergence angle [deg]    
        self.phi_div = phi_div * np.pi/180                          # divergence angle [deg]
        self.expansion_ratio = np.pi * ((D_e/2)**2 - (D_t/2)**2)    # expansion ratio (Ae/At)
        self.dx = step_size                                         # step size for the discretisation [-]

    def contour(self):
        # generate chamber countour x, y coordiantes
        
        # Cylindrical section
        x_cyl = np.arange(0, self.L_cyl, self.dx)
        y_cyl = np.ones(len(x_cyl)) * self.D_c / 2

        # Converging section radius
        n1 = round((self.r_1 * self.phi_conv ) / self.dx)     # number of discretisations in the radius
        t = np.linspace(0, self.phi_conv, n1)

        x_B = self.L_cyl + self.r_1 * np.sin(t)
        y_B = (self.D_c / 2 - self.r_1) + self.r_1 * np.cos(t)

        L_r1 = self.r_1 * np.sin(self.phi_conv)

        # Converging section 
        a_C = np.tan(self.phi_conv)
        L_con = ((self.D_c / 2 - self.r_1 + self.r_1 * np.cos(self.phi_conv)) - (self.D_t / 2 + self.r_2 - self.r_2 * np.cos(-self.phi_conv))) / a_C
        x_C = np.arange((self.L_cyl + L_r1 + self.dx), (self.L_cyl + L_r1 + L_con), self.dx)
        y_C = (self.D_c / 2 - self.r_1 + self.r_1 * np.cos(self.phi_conv)) - a_C * np.arange(self.dx, L_con, self.dx);


        # Converging radius up to throat
        n2 = round((self.r_2 * self.phi_conv) / self.dx)      # number of discretisations in the radius
        t = np.linspace(- self.phi_conv, 0, n2)

        L_r21 = self.r_2 * np.sin(self.phi_conv)
        x_D = self.L_cyl + L_r1 + L_con + L_r21 + self.r_2 * np.sin(t)
        y_D = self.D_t / 2 + self.r_2 - self.r_2 * np.cos(t)


        # Throat
        n3  = round((self.r_2 * self.phi_div) / self.dx)       # number of discretisations in the radius
        t = np.linspace(self.phi_div / n3, self.phi_div, n3)

        L_r22 = self.r_2 * np.sin(self.phi_div)

        x_E = self.L_cyl + L_r1 + L_con + L_r21 + self.r_2 * np.sin(t)
        y_E = self.D_t / 2 + self.r_2 - self.r_2 * np.cos(t)

        # Conical diverging secion 
        a_F = np.tan(self.phi_div)
        L_div=((self.D_e / 2) - (self.D_t / 2 + self.r_2 - self.r_2 * np.cos(self.phi_div))) / a_F
        x_F = np.arange((self.L_cyl + L_r1 + L_con + L_r21 + L_r22 + self.dx), (self.L_cyl + L_r1 + L_con + L_r21 + L_r22 + L_div + self.dx), self.dx)
        y_F = self.D_t / 2 + self.r_2 - self.r_2 * np.cos(self.phi_div) + a_F * np.arange(self.dx, L_div + self.dx, self.dx)


        self.x = np.concatenate((x_cyl, x_B, x_C, x_D, x_E, x_F), axis=0)
        self.y = np.concatenate((y_cyl, y_B, y_C, y_D, y_E, y_F), axis=0)
 

        xy = [[self.x[i], self.y[i]] for i in range(len(self.x))]
        self.geometry = np.vstack(xy)


    def plot_geometry(self):
        plt.plot(self.x, self.y, color='b')
        plt.plot(self.x, -self.y, color='b')
        plt.xlabel('x_coordinate [m]')
        plt.ylabel('y_coordinate [m]')
        plt.axis('equal')
        plt.show()



class CoolingGeometry():
    def __init__(self, chamber_geometry, h_c, phi, t_w_i, t_w_o, n_channels):
        self.geom = chamber_geometry                        # points along the chamber contour [x,y]
        self.h_c = h_c                                      # height of the cooling channels [m]
        self.phi = phi                                      # ratio of coverage of the cooling channels to solid wall
        self.t_w_i = t_w_i                                  # inner chamber wall thickness [m] 
        self.t_w_o = t_w_o                                  # outer chamber wall thickness [m] 
        self.n_channels = n_channels                        # number of cooling channels
        self.psi_c = 2 * np.pi * phi / n_channels           # radial section of single cooling channel 
        self.psi_w = 2 * np.pi * (1 - phi) / n_channels     # raidal section of single uncooled wall segment  

    def channel_geometry(self):
        # generate the cooling channel geometry; area, hydraulic diameter and wall thickness
        self.r_i = self.geom[:,1] + self.t_w_i                                      # inner radius of cooling channels 
        self.r_o = self.geom[:,1] + self.t_w_i + self.h_c                           # outer radius of cooling channels 

        self.A_c = self.psi_c / (self.psi_c + self.psi_w) * np.pi * (self.r_o**2 - self.r_i**2)     
        self.A_fin = self.psi_w / (self.psi_c + self.psi_w) * np.pi * (self.r_o**2 - self.r_i**2)                                              # cooling channel area [m^2]
        self.U_c = self.psi_c / (self.psi_c + self.psi_w) * 2 * np.pi * (self.r_o + self.r_i) + 2 * self.n_channels * (self.r_o + self.r_i)      # cooling channel circumference [m]
        self.D_h = 4 * self.A_c / self.U_c                                           # cooling channel hydraulic diameter [m]
        
        self.t_w_i_arr = np.ones(len(self.r_i)) * self.t_w_i                         # array of constant wall thickness
        self.t_w_o_arr = np.ones(len(self.r_i)) * self.t_w_o  

        self.U_cc = self.geom[:,1] * (self.psi_c + self.psi_w)                       # combustion chamber section circumference 
        self.U_c_i = self.r_i  * self.psi_c                                          # cooling channel innner side circumference
        self.U_c_s = self.r_o - self.r_i                                             # cooling channel side circumference
        self.U_c_o = self.r_o  * self.psi_c                                          # cooling channel outer side circumference

    def channel_efficiency(self, k, halpha_c, idx):
        # estimate of cooling channel efficiency based on heat trasnfer coefficient and channel shape 
        # k :                   [float] thermal conductivity 
        # halpha_c :            [float] convective heat transfer coefficient of the cooling channel
        # idx :                 [int]   index denoting location along chamber contour
        # -> halpha_c_corr:     [float] corrected convective heat transfer coefficient of the cooling channel

        a_c = self.r_i[idx] * self.psi_c
        b_c = self.r_i[idx] * self.psi_w

        eta_c = np.tanh(np.sqrt((2 * halpha_c * b_c ) / k) * self.t_w_i / b_c) / (np.sqrt((2 * halpha_c * b_c ) / k) * self.t_w_i / b_c)
  
        halpha_c_corr = (a_c + eta_c * (2 * (self.r_o[idx] - self.r_i[idx]) + a_c)) / (a_c + b_c) * halpha_c
    
        return halpha_c_corr        



if __name__ == "__main__":
    # test functionality of the Chamber and Cooling Geometry Classes

    D_c = 0.014
    D_t = 0.0042
    D_e = 0.0052
    L_cyl = 0.05
    r_1 = 0.01203
    r_2 = 0.003
    phi_conv = 30
    phi_div = 15

    cg = ChamberGeometry(D_c, D_t, D_e, L_cyl, r_1, r_2, phi_conv, phi_div)
    cg.contour()
    cg.plot_geometry()

    n = 8 
    h_c = 3e-3
    phi = 1/3
    t_w_i = 1e-3
    t_w_o = 1e-3

    c = CoolingGeometry(cg.geometry, h_c, phi, t_w_i, t_w_o, n)
    c.channel_geometry()


