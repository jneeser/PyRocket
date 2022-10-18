import numpy as np
import config

# output class for 1D and 2D sims

class Output1D():
    def __init__(self):
        self.halpha         = np.ndarray(len(config.geometry[:,1]))
        self.halpha_c       = np.ndarray(len(config.geometry[:,1]))
        self.q_rad          = np.ndarray(len(config.geometry[:,1]))
        self.q              = np.ndarray(len(config.geometry[:,1]))
        self.T_wall_i       = np.ndarray(len(config.geometry[:,1]))
        self.T_wall_o       = np.ndarray(len(config.geometry[:,1]))
        self.T_c            = np.ndarray(len(config.geometry[:,1]))
        self.P_c            = np.ndarray(len(config.geometry[:,1]))
        self.Re             = np.ndarray(len(config.geometry[:,1]))
        self.v_coolant      = np.ndarray(len(config.geometry[:,1]))
        self.T_hg           = np.ndarray(len(config.geometry[:,1]))