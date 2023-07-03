# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

import numpy as np
import matplotlib.pyplot as plt
from ActiveMaterial import ActiveMaterial
from V_waveform import V_waveform
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
            dst = r'C:\Users\aldanads\OneDrive - Trinity College Dublin\2D device simulator project\Publications\Schottky barrier\Simulations\Testing\\'
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
         Davelou, Daphne, Georgios Kopidakis, George Kioseoglou, and Ioannis N. Remediakis. 
         "MoS2 nanostructures: Semiconductors with metallic edges." 
         Solid state communications 192 (2014): 42-46.
         er = 7.1 for bulk --> DFT calculations
         
         Santos, Elton JG, and Efthimios Kaxiras. 
         "Electrically driven tuning of the dielectric constant in MoS2 layers." 
         ACS nano 7, no. 12 (2013): 10741-10746.
         er = 4-16 for 8L under electric field
         
         Paul, Atanu, and Ilya Grinberg. 
         "Optical Properties of Mo S 2, Mo Se 2, WS 2, and W Se 2 under External Electric Field." 
         Physical Review Applied 17, no. 2 (2022): 024042.
         er = 11.5 for 5L under electric field

    """
    er = 21 # Relative permitivity
    screening = 0.05 # Screening of electric field
    q = -2 # Particle charge
    
    """
     Electrodes properties
    """
    # Position where finish or start the electrode
    electrodes = (round(device_size[0]/(steps[0]*6)),round(5*device_size[0]/(steps[0]*6))) 

    
    """
    # =============================================================================
    #     # Activation energies 
    # =============================================================================
    """
    E_diffusion = 1.1
    
    # Defects can migrate up-down, lateral x-axis and y-axis
    Act_energy = [E_diffusion,E_diffusion,E_diffusion,E_diffusion,E_diffusion,E_diffusion]
    
    # Temperature in kelvin (k)
    T = 300
    
    """
     Sangwan, Vinod K., Hong-Sub Lee, Hadallia Bergeron, Itamar Balla, Megan E. Beck, Kan-Sheng Chen, and Mark C. Hersam. 
     "Multi-terminal memtransistors from polycrystalline monolayer molybdenum disulfide." Nature 554, no. 7693 (2018): 500-504.
     Sangwan says that phi_bn (effective Schottky barrier) range between 80–125 meV in experiments
     They use 20 meV to 280 meV in simulations
     """
    phi_b0 = 0.385 # Schottky barrier (eV) not modulated  (fitting parameter)
    w = 3E-9 # w (nm) is in the range of a few nanometers --> region with excess of dopants (fitting parameter)
    A = 10E-5 # A constant
    R = 5000
    ideality_factor = 10
    tol = 1e-8
    
    square_size = 4
    # Introducing defects in the active layer
    n_defects = 50
    n = n_defects/np.product(device_size) # Density of defects in the structure

    
    schottky_parameters = [phi_b0, w, A, n, R, ideality_factor,tol]
    
    
    # Voltage waveform
    shape_phase = ['forming','triangle']
    MaxV = -1 
    wave_width = 4
    delta_t = 0.1 # Sampling time
    V_forming = 3
    t_forming = 3
    V = V_waveform(shape_phase,MaxV,wave_width,delta_t,V_forming,t_forming)
    
    MoS2_layer = ActiveMaterial(Act_energy,device_size,steps,time,er,screening,electrodes,q,T,schottky_parameters,square_size)


    defects_list = MoS2_layer.introducing_defects(rng, n_defects)
    #defects_list = MoS2_layer.test_1(rng)
    
    
    return MoS2_layer,paths,rng,defects_list,V

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
            
        