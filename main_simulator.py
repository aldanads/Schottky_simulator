# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:52:56 2023

@author: ALDANADS
"""
from initialization import *
from KMC import KMC
import time
import matplotlib.pyplot as plt


end = time.time()

save_data = False
save_var = False

n_sim = 1
MoS2_layer,paths,rng,defects_list,V = initialization(n_sim,save_data)

carry_on = True
i = 0

current = []
elapsed_time = []

MoS2_layer.plot_particles()
V.update_V(MoS2_layer.time[-1])
start = time.time()
MoS2_layer.SolvePotentialAndField(V.voltage[-1])
end = time.time()
elapsed_time.append(end-start)

current.append(MoS2_layer.Schottky_current(V.voltage[-1]))



while carry_on:
    i += 1
    V.update_tmax(i)
    start = time.time()
    
    while MoS2_layer.time[-1] < V.tmax:
        MoS2_layer,defects_list = KMC(MoS2_layer,rng,defects_list,V)
        
    MoS2_layer.plot_particles()
    V.update_V(MoS2_layer.time[-1])
    MoS2_layer.SolvePotentialAndField(V.voltage[-1])
    current.append(MoS2_layer.Schottky_current(V.voltage[-1]))

    print(f'Voltage (V): {V.voltage[-1]:.2f}', f'Time (s): {MoS2_layer.time[-1]:.2f}',f'Current (A): {current[-1]:.4e}')
    plt.semilogy(V.voltage,abs(np.array(current)))
    plt.ylim(1e-11,1e-3)
    end = time.time()
    elapsed_time.append(end-start)

    if i > 1000:
        break
    