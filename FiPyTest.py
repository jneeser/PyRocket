import config
import RegenCooling as rc


from fipy import CellVariable, Gmsh2D, DiffusionTerm, Viewer, Grid2D, TransientTerm
from fipy.tools import numerix
from fipy import input
from scipy import interpolate
import numpy as np


def get_mesh(cellSize):
    x = config.geometry[-1,0]
    y = config.t_w_i

    mesh = Gmsh2D('''
              cellSize = %(cellSize)g;
              x = %(x)g;
              y = %(y)g;
              Point(1) = {0, 0, 0, cellSize};
              Point(2) = {x, 0, 0, cellSize};
              Point(3) = {x, y, 0, cellSize};
              Point(4) = {0, y, 0, cellSize};
              Line(6) = {1, 2};
              Line(7) = {2, 3};
              Line(8) = {3, 4};
              Line(9) = {4, 1};
              Line Loop(10) = {6, 7, 8, 9};
              Plane Surface(11) = {10};
              ''' % locals()) 

    return mesh

cellSize = 0.0008
mesh = get_mesh(cellSize)
phi = CellVariable(name = "T [K]",
                   mesh = mesh,
                   value = 500.)

if __name__ == '__main__':
    viewer = Viewer(vars=phi, datamin=288, datamax=1400)

    viewer.plotMesh()
    x, y = mesh.faceCenters

    method = 'cinjarew'
    cooling_method = 'gnielinski'
    sim = rc.HeatTransfer(config.cea, config.gas, config.chamber.geometry, config.material, config.coolant, config.cooling_geom, config.m_dot, config.m_dot_coolant, method, cooling_method, config.correction, config.eta_combustion)
    sim.run_sim()

    def halpha_gas(x):
        halpha = interpolate.interp1d(config.geometry[:,0], sim.out.halpha, kind='linear', fill_value='extrapolate')
        return halpha(x)

    def halpha_coolant(x):
        halpha_c = interpolate.interp1d(config.geometry[:,0], sim.out.halpha_c,                                                         kind='linear', fill_value='extrapolate')
        return halpha_c(x)

    def coolant_temp(x):
        T_c = interpolate.interp1d(config.geometry[:,0], sim.out.T_c, kind='linear', fill_value='extrapolate')
        return T_c(x)

    def hot_gas_temp(x):
        T_aw = interpolate.interp1d(config.geometry[:,0], sim.gas.T_aw, kind='linear', fill_value='extrapolate')
        return T_aw(x)

    def thermal_conductivity(T):
        k = config.material.k
        return k(T)



    D = 1e-4

    topFlux     = halpha_coolant(x) * mesh.faceNormals
    bottomFlux  = halpha_gas(x) * mesh.faceNormals
    coolantTemp = coolant_temp(x) * mesh.faceNormals
    hotGasTemp  = - hot_gas_temp(x) * mesh.faceNormals
    k           = thermal_conductivity(phi.faceValue)


    phi.faceGrad.constrain([-1/k * topFlux * (phi.faceValue - coolantTemp)], mesh.facesTop)
    phi.faceGrad.constrain([-1/k * bottomFlux * (phi.faceValue - hotGasTemp)], mesh.facesBottom)

    eq = TransientTerm() == DiffusionTerm(coeff=D)

    timeStepDuration = 10 * 0.99 * (cellSize)**2 / (2*D)

    steps = 200
    t = 0

    print('Time [s]', '     ', 'max T [K]')

    for step in range(steps):

        eq.solve(var=phi, dt=timeStepDuration)
        t += timeStepDuration

        print(round(t,3), '   ', round(max(phi),3))

        viewer.plot()

