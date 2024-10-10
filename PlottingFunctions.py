#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3			                #	 
#                                                                   #
#####################################################################

import matplotlib.pyplot as plt
import numpy as np

import config

class Plotting1D():
    def __init__(self, save_path, save=True, show=False):
        self.save_path = save_path
        self.save = save
        self.show = show
        self.data = np.genfromtxt(save_path+'/sim_data.csv', delimiter=",", skip_header=1)


    def temperature_plot(self):
        plt.plot(config.geometry[:,0]*1e3, self.data[:,7]-273.15, color='b', label='Coolant bulk temperature')
        plt.plot(config.geometry[:,0]*1e3, self.data[:,6]-273.15, color='r', linestyle='--', label='Inner wall temperature')
        plt.xlabel('x coordiante [mm]')
        plt.ylabel('Temperature [$^\circ$C]')
        plt.legend(loc='best')

        if self.save:
            plt.savefig(self.save_path + '/temperature', dpi=600)
            plt.close()
            
        if self.show:
            plt.show()
        

    def pressure_plot(self):
        plt.plot(config.geometry[:,0]*1e3, self.data[:,8]/1e5, color='b', label='Coolant pressure')
        plt.xlabel('x coordiante [mm]')
        plt.ylabel('Pressure [bar]')
        plt.legend(loc='best')

        if self.save:
            plt.savefig(self.save_path + '/pressure', dpi=600)
            plt.close()
        
        if self.show:
            plt.show()
        

    def heat_transfer_coeff_plot(self):
        plt.plot(config.geometry[:,0]*1e3, self.data[:,2], color='b', label='$h_{alpha}$ gas')
        plt.plot(config.geometry[:,0]*1e3, self.data[:,3], color='r', linestyle='--',  label='$h_{alpha}$ coolant')
        plt.xlabel('x coordiante [mm]')
        plt.ylabel('Heat transfer coefficient [W/m^2/K]')
        plt.legend(loc='best')
        
        if self.save:
            plt.savefig(self.save_path + '/halpha', dpi=600)
            plt.close()

        if self.show:
            plt.show()
    

    def heat_flux_plot(self):
        plt.plot(config.geometry[:,0]*1e3, self.data[:,4]/1e6, color='b', label='q radiation')
        plt.plot(config.geometry[:,0]*1e3, self.data[:,5]/1e6, color='r', linestyle='--', label='q total')
        plt.xlabel('x coordiante [mm]')
        plt.ylabel('Heat flux [MW/m^2]')
        plt.legend(loc='best')
        
        if self.save:
            plt.savefig(self.save_path + '/heat_flux', dpi=600)
            plt.close()
        
        if self.show:
            plt.show()
    

    def reynolds_plot(self):
        plt.plot(config.geometry[:,0]*1e3, self.data[:,9], color='b', label='Coolant Re')
        plt.xlabel('x coordiante [mm]')
        plt.ylabel('Re')
        plt.legend(loc='best')

        if self.save:
            plt.savefig(self.save_path + '/reynolds', dpi=600)
            plt.close()
    
        if self.show:
            plt.show()
        

def multi_plot(data1, data2, data3, data4, label1, label2, label3, label4):
    f, axes = plt.subplots(4, 1)
    #f.subplots_adjust(hspace = 0)

    axes[0].plot(config.geometry[:,0], data1)
    axes[0].set_ylabel(label1)

    axes[1].plot(config.geometry[:,0], data2)
    axes[1].set_ylabel(label2)

    axes[2].plot(config.geometry[:,0], data3)
    axes[2].set_ylabel(label3)

    axes[3].plot(config.geometry[:,0], data4)
    axes[3].set_ylabel(label4)

    plt.xlabel('x coordiante [m]')
    plt.show()


if __name__ == "__main__":
    plot = Plotting1D(save_path=config.save_path, save=True, show=True)
    plot.temperature_plot()
    plot.pressure_plot()
    plot.heat_transfer_coeff_plot()
    plot.heat_flux_plot()
    plot.reynolds_plot()
