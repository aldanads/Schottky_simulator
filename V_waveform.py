# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:28:01 2023

@author: ALDANADS
"""
import numpy as np


class V_waveform():
    
    def __init__(self,shape_phase,MaxV,wave_width,delta_t,V_forming,t_forming,n_cycles):
        
        self.shape = shape_phase[0]
        self.shape_phase = shape_phase
        self.MaxV = MaxV
        self.V_forming = V_forming
        self.t_forming = t_forming
        self.wave_width = wave_width # waveform width (s)
        self.delta_t = delta_t # sampling time (s)
        self.voltage = []
        self.tmax = 0
        self.cycles = 0
        self.n_cycles = n_cycles
        
    def update_V(self,t):
        
        if self.shape == 'forming':
            self.forming(t)
            
        elif self.shape == 'triangle':
            self.triangle(t)
        
    def triangle(self,t):
    
        MaxV = self.MaxV
        wave_width = self.wave_width
        dVdt = (2*MaxV)/wave_width # Rate of change of voltage
        
        T = wave_width * 2 # Period of the voltage waveform
            
        phase = (t-self.t_forming-self.delta_t) % T
        quadrant = np.floor((phase*4) / T)
        q_phase = phase - (quadrant * (wave_width/2))
        
        if quadrant == 0:
            V = q_phase * dVdt
        elif quadrant == 1:
            V = MaxV - (dVdt * q_phase)
        elif quadrant == 2:
            V = -dVdt * q_phase
        elif quadrant == 3:
            V = - MaxV + (dVdt * q_phase)
            
        if phase == 0:
            self.cycles += 1
            
        self.voltage.append(V)
        
    def forming(self,t):
        
        V_forming = self.V_forming
        t_forming = self.t_forming
        
        dVdt = V_forming/t_forming
                        
        if t >= t_forming:
            self.shape = self.shape_phase[1]
        
        self.voltage.append(t * dVdt)

        
    def update_tmax(self,i):
        
        self.tmax = i * self.delta_t

            
        
        

    