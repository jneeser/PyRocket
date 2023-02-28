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


# aniline
card_str = """
fuel C6H7N(L)  C 6.0   H 7.0    N 1.0  wt%=100
h, kj/mol=31.3    t(k)=298.15   rho=1.03 
"""
add_new_fuel( 'aniline', card_str )
aniline = rocketcea.blends.newFuelBlend(fuelL=['aniline'], fuelPcentL=[100]) 


# furfurylalcohol
card_str = """
fuel C5H6O2(L)  C 5.0   H 6.0    O 2.0  wt%=100
h, kj/mol=-276.2    t(k)=298.15   rho=1.13 
"""
add_new_fuel( 'furfurylalcohol', card_str )
furfurylalcohol = rocketcea.blends.newFuelBlend(fuelL=['furfurylalcohol'], fuelPcentL=[100]) 



