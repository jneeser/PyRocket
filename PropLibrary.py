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
import rocketcea
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant



# Hydrogen Peroxide  
peroxide98 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[98,2]) 
peroxide96 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[96,4]) 
peroxide94 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[94,6]) 
peroxide85 = rocketcea.blends.newOxBlend(oxL=['H2O2', 'H2O'], oxPcentL=[85,15]) 


# EMIM SCN 
card_str = """
fuel C7H11N3S(L)  C 7.0   H 11.0   N 3.0   S 1.0  wt%=100
h, kj/mol=52.8    t(k)=298.15   rho=1.11 
"""
add_new_fuel( 'EMIMSCN', card_str )
EMIMSCN = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN'], fuelPcentL=[100]) 

