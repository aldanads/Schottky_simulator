# -*- coding: utf-8 -*-
"""
Created on Mon May 22 19:46:20 2023

@author: ALDANADS
"""
import numpy as np
from balanced_tree import Node, build_tree, update_data, search_value


def KMC(MoS2_layer,rng,defects_list):
    
    
    TR_list = [] # Catalog of  non-zero transition rates in the system
    
       
    # Collect non-zero transition rates for each defect
    for i, defect in enumerate(defects_list):
        # Transition rate of event j of defect i
        TR_list.extend([(TR, j, i) for j, TR in enumerate(defect.TR) if TR != 0.0])
    

    chosen_event,time = choose_event(TR_list,rng) # Choose an event and time step
    
    Vs = defects_list[chosen_event[2]]  # Get the defect corresponding to the chosen event
    Vs.processes(chosen_event[1], MoS2_layer)  # Apply the chosen event to the system

    MoS2_layer.time_event_track(time, chosen_event[1])  # Update time and event count
    
    
    return MoS2_layer,defects_list



"""
 ---- chosen_event,time = choose_event(TR_list,rng) ------
 
# Organize transition rates in a balanced binary tree
# TR_list is a list of tuples(TR,event_j,particle_i)
# Return the chosen event and the time step corresponding to this kMC step
"""
def choose_event(TR_list,rng):
    
    # Sort the list of events
    TR_list.sort(key=lambda x: x[0])
    
    TR_tree = build_tree(TR_list) # Build a balanced tree structure
    
    # Each node is the sum of their children, starting from the leaf
    sumTR = update_data(TR_tree) # Head node is sum of all transition rates

    # When there is only one node it returns a tuple
    sumTR = sumTR[0] if isinstance(sumTR, tuple) else sumTR
    
    # Remember tuple (TR,j_event,i_particle)
    chosen_event = search_value(TR_tree,sumTR*rng.random()) # Find the chosen event

    #Calculate the time step
    time = -np.log(rng.random()) / sumTR
    
    return chosen_event,time
    