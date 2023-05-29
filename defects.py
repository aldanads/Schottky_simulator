# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:56:36 2023

@author: ALDANADS
"""

import numpy as np

class Defects():
    
    
    def __init__(self,coord_xyz,Act_energy,Grid_states):
        
        self.coord_xyz = coord_xyz
        self.Act_energy = Act_energy
        self.length_x = len(Grid_states)-1
        self.length_y = len(Grid_states[0])-1
        self.length_z = len(Grid_states[0][0])-1
        
        self.neighbors(Grid_states)

        self.transition_rates()
        
    
    """
    ---------------- Calculate neighbors -------------------------
    """   
    def neighbors(self,Grid_states):
        
        List_neighbors = [0] * 6 # Six directions: 2 each axix
 
        i,j,k = self.coord_xyz[0],self.coord_xyz[1],self.coord_xyz[2]

        
        # Checking if there is a particle or a domain boundary
        # Check lateral x-axis
        if (i == self.length_x) or (Grid_states[i+1,j,k] == 1):
            List_neighbors[0] = (i+1,j,k)
            
            
        if (i == 0) or (Grid_states[i-1,j,k] == 1):
            List_neighbors[1] = (i-1,j,k)
            
        # Check lateral y-axis
        if (j == self.length_y) or (Grid_states[i,j+1,k] == 1):
            List_neighbors[2] = (i,j+1,k)
            
        if (j == 0) or (Grid_states[i,j-1,k] == 1):
            List_neighbors[3] = (i,j-1,k)
            
        #Check up and down
        if (k == self.length_z) or (Grid_states[i,j,k+1] == 1):
            List_neighbors[4] = (i,j,k+1)
            
        if (k == 0) or (Grid_states[i,j,k-1] == 1):
            List_neighbors[5] = (i,j,k-1)
            
        self.List_neighbors = List_neighbors
        
    """
    ---------------- Calculate Transition Rates -----------------
    """ 
    def transition_rates(self,Ex=None,Ey=None,Ez=None,steps=None):
        
        i,j,k = self.coord_xyz[0],self.coord_xyz[1],self.coord_xyz[2]
        
        energy_modulation = 0
        if Ex is not None:
            # Energy modulation by the electric field - These indexes match Act_energy indexes
            energy_modulation = np.array([Ex[i,j,k]*steps[0],-Ex[i,j,k]*steps[0],
                                 Ey[i,j,k]*steps[1],-Ey[i,j,k]*steps[1],
                                 Ez[i,j,k]*steps[2],-Ez[i,j,k]*steps[1]])*1e-9
            
        
        available_idx = [i for i, item in enumerate(self.List_neighbors) if not isinstance(item, tuple)]
        TR = np.zeros(len(self.Act_energy))
        
        kb = 8.6173324E-5 # Boltzmann constant
        nu0=7E13;  # nu0 (s^-1) bond vibration frequency
        T = 300
        Energy_mod = (np.array(self.Act_energy)-energy_modulation)
        Energy_mod[Energy_mod < 0.1] = 0.1

        TR[available_idx] = nu0*np.exp(-Energy_mod[available_idx]/(kb*T))
      
        
        self.TR = TR
    

        

    
    """
    ---------------- Apply event -----------------
    """ 
    def processes(self,selected_event,MoS2_layer):
        
        i,j,k = self.coord_xyz[0],self.coord_xyz[1],self.coord_xyz[2]
        
        # Diffusion in lateral x-axis
        if selected_event == 0:
            MoS2_layer.Grid_states[i,j,k] = 0
            MoS2_layer.Grid_states[i+1,j,k] = 1
            i += 1
            
        elif selected_event == 1:
            MoS2_layer.Grid_states[i,j,k] = 0
            MoS2_layer.Grid_states[i-1,j,k] = 1
            i -= 1
            
        # Diffusion in lateral y-axis
        elif selected_event == 2:
            MoS2_layer.Grid_states[i,j,k] = 0
            MoS2_layer.Grid_states[i,j+1,k] = 1
            j += 1
            
        elif selected_event == 3:
            MoS2_layer.Grid_states[i,j,k] = 0
            MoS2_layer.Grid_states[i,j-1,k] = 1
            j -= 1
            
        # Diffusion in lateral y-axis
        elif selected_event == 4:
            MoS2_layer.Grid_states[i,j,k] = 0
            MoS2_layer.Grid_states[i,j,k+1] = 1
            k += 1
            
        elif selected_event == 5:
            MoS2_layer.Grid_states[i,j,k] = 0
            MoS2_layer.Grid_states[i,j,k-1] = 1
            k -= 1
            
        self.coord_xyz[0],self.coord_xyz[1],self.coord_xyz[2] = i,j,k
        self.neighbors(MoS2_layer.Grid_states)
        self.transition_rates()
    
            
