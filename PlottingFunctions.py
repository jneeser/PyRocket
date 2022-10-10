import matplotlib.pyplot as plt
import numpy as np
import scipy

import config


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

