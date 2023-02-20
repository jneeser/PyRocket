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

# EMIM SCN 
card_str = """
fuel C7H11N3S(L)  C 7.0   H 11.0   N 3.0   S 1.0  wt%=100
h, kj/mol=52.8    t(k)=298.15   rho=1.11 
"""
add_new_fuel( 'EMIMSCN', card_str )
EMIMSCN = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN'], fuelPcentL=[100]) 

# CuSCN 
card_str = """
fuel C2CuN2S2(L)  C 2.0   N 3.0   S 2.0  Cu 1.0  wt%=100
h, kj/mol=-97.97    t(k)=298.15   rho=2.88
"""
add_new_fuel( 'CuSCN', card_str )
CuSCN = rocketcea.blends.newFuelBlend(fuelL=['CuSCN'], fuelPcentL=[100]) 

# HimSCN 
card_str = """
fuel C4H5N3S(L)  C 4.0   H 5.0   N 3.0   S 1.0  wt%=100
h, kj/mol=127.0    t(k)=298.15   rho=1.355 
"""
add_new_fuel( 'HimSCN', card_str )
HimSCN = rocketcea.blends.newFuelBlend(fuelL=['HimSCN'], fuelPcentL=[100]) 

# Ionic Fuel Blends
EMIMSCN = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN'], fuelPcentL=[100]) 
HIP11 = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN', 'CuSCN'], fuelPcentL=[95,5]) 
HIM35 = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN', 'HimSCN'], fuelPcentL=[65,35]) 
HIM30_Cu5 = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN', 'HimSCN',  'CuSCN'], fuelPcentL=[65,30,5]) 
HIM25_Cu5 = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN', 'HimSCN',  'CuSCN'], fuelPcentL=[70,25,5]) 
HIM30_Cu3 = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN', 'HimSCN',  'CuSCN'], fuelPcentL=[67,30,3]) 
HIM25_Cu3 = rocketcea.blends.newFuelBlend(fuelL=['EMIMSCN', 'HimSCN',  'CuSCN'], fuelPcentL=[72,25,3]) 

