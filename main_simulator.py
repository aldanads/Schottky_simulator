# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:52:56 2023

@author: ALDANADS
"""
from initialization import *
from KMC import KMC

save_data = False
save_var = True

n_sim = 1
MoS2_layer,paths,rng,defects_list = initialization(n_sim,save_data)

carry_on = True
i = 0
while carry_on:
    i += 1
    
    MoS2_layer,defects_list = KMC(MoS2_layer,rng,defects_list)
    #MoS2_layer.SolvePotentialAndField(2)
    MoS2_layer.plot_particles()
    
    if i > 50:
        break
    