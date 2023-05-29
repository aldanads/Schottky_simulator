# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

import numpy as np
import matplotlib.pyplot as plt
from ActiveMaterial import ActiveMaterial
import shutil
import os 
import platform



def initialization(n_sim,save_data):

    # Random seed as time
    rng = np.random.default_rng() # Random Number Generator (RNG) object

    # Default resolution for figures
    plt.rcParams["figure.dpi"] = 300 # Default value of dpi = 300
    
    if save_data:
        files_copy = ['defects.py', 'hex_lattice.py', 'initialization.py','KMC.py','main_simulator.py','load_variables.py']
        
        if platform.system() == 'Windows': # When running in laptop
            dst = r'C:\Users\aldanads\OneDrive - TCDUD.onmicrosoft.com\2D device simulator project\Publications\Layer growth\Simulations\kMC2\Horizontal\\'
        elif platform.system() == 'Linux': # HPC works on Linux
            dst = r'/home/users/aldanads/Crystal growth/Simulations/Growth rate -  Diffusion rate/1.2eV/'
            
        paths = save_simulation(files_copy,dst,n_sim) # Create folders and python files
    else:
        paths = {'data': ''}
        
    

    """
     Grid parameters
    """
    steps = (0.639,0.639,0.639) # nm | x-y-z
    device_size = (20, 20, 20) # Size of the grid in nm. grid_size[0] is x, grid_size[1] is y and grid_size[3] is z.
    
    time = 0 # Initialization of time
    
    """
     Electrical properties
    """
    er = 21 # Relative permitivity
    screening = 0.025 # Screening of electric field
    q = 2 # Particle charge
    
    """
     Electrodes properties
    """
    # Position where finish or start the electrode
    electrodes = (round(device_size[0]/(steps[0]*8)),round(7*device_size[0]/(steps[0]*8))) 

    
    """
    # =============================================================================
    #     # Activation energies 
    # =============================================================================
    """
    E_diffusion = 1.5
    
    # Defects can migrate up-down, lateral x-axis and y-axis
    Act_energy = [E_diffusion,E_diffusion,E_diffusion,E_diffusion,E_diffusion+0.5,E_diffusion+0.5]
    
    # Temperature in kelvin
    T = 300
    

    
    
    
    
    
    
    MoS2_layer = ActiveMaterial(Act_energy,device_size,steps,time,er,screening,electrodes,q)

    # Introducing defects in the active layer
    n_defects = 50
    defects_list = MoS2_layer.introducing_defects(rng, n_defects)
    #defects_list = MoS2_layer.test_1(rng)
    
    
    return MoS2_layer,paths,rng,defects_list

def save_simulation(files_copy,dst,n_sim):
    

    if platform.system() == 'Windows':
        parent_dir = 'Sim_'+str(n_sim)+'\\'
        os.makedirs(dst+parent_dir) 
        dst = dst+parent_dir
        program_directory = 'Program\\'
        data_directoy = 'Crystal growth\\'
        
    elif platform.system() == 'Linux':
        parent_dir = 'Sim_'+str(n_sim)+'/'
        os.makedirs(dst+parent_dir) 
        dst = dst+parent_dir
        program_directory = 'Program/'
        data_directoy = 'Crystal growth/'

    os.makedirs(dst + program_directory)
    os.makedirs(dst + data_directoy)
    
    paths = {'data': dst + data_directoy, 'program': dst + program_directory}

    for files in files_copy:
        shutil.copyfile(files, paths['program']+files)
        
    return paths

def save_variables(paths,variables):
    
    
    if platform.system() == 'Windows': # When running in laptop

        import shelve
    
        filename = 'variables'
        my_shelf = shelve.open(paths+filename,'n') # 'n' for new
        
        for key in variables:
            my_shelf[key] = variables[key]
    
        my_shelf.close()

    elif platform.system() == 'Linux': # HPC works on Linux
    
        import pickle
    
        filename = 'variables.pkl'    
    
        # Open a file and use dump()
        with open(paths+filename, 'wb') as file:
              
            # A new file will be created
            pickle.dump(variables,file)
            
        