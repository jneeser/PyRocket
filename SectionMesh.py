#####################################################################
#                           PyRocket								#
# 2D Regenertive Cooling Simulation for Bipropellant Rocket Engines #
#                                                                   #
# Creator:  Joanthan Neeser                                         #
# Date:     15.12.2022                                              #
# Version:  2.3  													#
# License:	GNU GENERAL PUBLIC LICENSE V3							#                                          
#####################################################################

import config
from fipy import Gmsh2D


def get_section_mesh(cell_size, idx):
    r_c = config.geometry[idx,1]            # inner chamber radius
    h_c = config.cooling_geom.h_c[idx]      # height of the cooling channel
    t_w_i = config.cooling_geom.t_w_i[idx]  # inner wall thickness  
    t_w_o = config.cooling_geom.t_w_o[idx]  # outer wall thickness
    psi_c = config.cooling_geom.psi_c[idx]       # radial section of Single cooling channel  
    psi_w = config.cooling_geom.psi_w[idx]       # radial section of Single wall segment 

    mesh = Gmsh2D('''
              cellSize = %(cell_size)g;
              r_c   = %(r_c)g;
              h_c   = %(h_c)g;
              t_w_i = %(t_w_i)g;
              t_w_o = %(t_w_o)g;
              psi_c = %(psi_c)g;
              psi_w = %(psi_w)g;

              Point(1) = {0, 0, 0};             
              
              Point(2) = {0, r_c, 0, cellSize};
              Point(9) = {0, r_c + t_w_i, 0, cellSize};
              Point(6) = {0, r_c + t_w_i + h_c, 0, cellSize};
              Point(5) = {0, r_c + t_w_i + h_c + t_w_o, 0, cellSize};
              
              Point(8) = {Sin(psi_c/2)*(r_c + t_w_i), Cos(psi_c/2)*(r_c + t_w_i), 0, cellSize};
              Point(7) = {Sin(psi_c/2)*(r_c + t_w_i + h_c), Cos(psi_c/2)*(r_c + t_w_i + h_c), 0, cellSize};
              
              Point(3) = {Sin(psi_c/2 + psi_w/2)*r_c, Cos(psi_c/2 + psi_w/2)*r_c, 0, cellSize};
              Point(4) = {Sin(psi_c/2 + psi_w/2)*(r_c + t_w_i + h_c + t_w_o), Cos(psi_c/2+ psi_w/2)*(r_c + t_w_i + h_c + t_w_o), 0, cellSize};
              
              Circle(11) = {2, 1, 3};
              Line(12)   = {3, 4};
              Circle(13) = {4, 1, 5};
              Line(14)   = {5, 6};
              Circle(15) = {6, 1, 7};
              Line(16)   = {7, 8};
              Circle(17) = {8, 1, 9};
              Line(18)   = {9, 2};

              Line Loop(20) = {11, 12, 13 , 14, 15, 16, 17, 18};
              Plane Surface(21) = {20};
              Physical Surface("Domain")         = {21};
              Physical Line("ChamberWall")       = {11};
              Physical Line("SideWall")          = {12};
              Physical Line("OuterWall")         = {13};
              Physical Line("TopWall")           = {14};
              Physical Line("CoolantTopWall")    = {15};
              Physical Line("CoolantSideWall")   = {16};
              Physical Line("CoolantBottomWall") = {17};
              Physical Line("BottomWall")        = {18};
              
              ''' % locals()) 

    return mesh
