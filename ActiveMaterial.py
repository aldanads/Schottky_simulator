# -*- coding: utf-8 -*-
"""
Created on Tue May 23 17:57:25 2023

@author: ALDANADS
"""
import numpy as np
from defects import Defects
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection




class ActiveMaterial():
    
    def __init__(self,Act_energy,device_size,steps,time,er,screening,electrodes,q,T,schottky_parameters,square_size):
        
        self.Act_energy = Act_energy
        self.device_size = device_size
        self.steps = steps
        self.time = []
        self.time.append(time)
        self.Grid_states = self.square_grid()
        self.events = [0]*len(Act_energy)
        self. T = T
        
        # Electric field parameters
        self.er = er
        self.q = q
        self.screening = screening
        self.electrodes = electrodes
        self.u = np.zeros(self.Grid_states.shape)
        
        # Defect density calculation
        self.square_size = square_size
        
        # Electrical parameters
        # Change to cm --> F/m = C / (V * m) = C / (100 * V * cm)
        e0 = 8.854187817E-12  # Permittivity of free space F/m = C / (V * m) = (C^2 / (N*cm^2))
        e_charge = 1.60217663E-19 # Electron charge (C)
        self.e_er = e_charge/(self.er*e0)
        self.A0 = 1.202E6 # Richardson constant (A*(m^-2)*(K^-2))
        kb=8.6173324E-5  # Boltzmann constant (eV/K)
        self.kb_T = kb * T
        
        self.phi_b0 = schottky_parameters[0]
        self.w = schottky_parameters[1]
        self.A = schottky_parameters[2]
        self.n = schottky_parameters[3]
        self.R = schottky_parameters[4]
        self.ideality_factor = schottky_parameters[5]
        self.tol = schottky_parameters[6]
        
        # Volume = (size x of the electrode + region of excess of dopants w) * (size y of electrode) * w in z direction
        self.vol_electrode = (electrodes[0]*steps[0] + schottky_parameters[1]*1e9) * device_size[1] * schottky_parameters[1]*1e9
        # Cross sectional area perpendicular to the current --> Area of the electrode (m**2)
        # The electrode are in nm scale
        self.cross_sectional_area =  device_size[1] * electrodes[0]*steps[0] *1e-18


        
    """  
    ---------------------------------------------
    --------------- Define grid -----------------
    ---------------------------------------------
    Sulfur vacancies = 1
    ---------------------------------------------
    """
    def square_grid(self):
        
        # Create MoS2 crystal - Amorphous: Square grid
        grid_points = [int(np.round(size / step)) for size, step in zip(self.device_size, self.steps)]

        return np.zeros(grid_points, dtype=int)
        
    
    def introducing_defects(self, rng, n_defects):
        
        grid_points = self.Grid_states.shape # Initialize grid size
        defects_list = [] # Initialize defects_list
        
        chosen_coordinates = set() # I want unique coordinates

        while len(defects_list) < n_defects:
            Vs_ijk = rng.integers(0, grid_points, size=3)

            if tuple(Vs_ijk) not in chosen_coordinates:
                chosen_coordinates.add(tuple(Vs_ijk))
                defects_list.append(Defects(Vs_ijk, self.Act_energy, self.Grid_states,self.q))
                self.Grid_states[tuple(Vs_ijk)] = 1

        return defects_list
    
    def test_1(self,rng):
        
        defects_list = [] # Initialize defects_list

        x_axis_paralallel = [np.array([6,6,6]),np.array([5,6,6])] 
        y_axis_paralallel = [np.array([6,6,6]),np.array([6,5,6])]
        z_axis_paralallel = [np.array([6,6,6]),np.array([6,6,5])]
        single_particle = [np.array([6,6,6])]
        
        test_defects = [x_axis_paralallel,y_axis_paralallel,z_axis_paralallel,single_particle]
        
        for Vs in test_defects[3]:
            defects_list.append(Defects(Vs, self.Act_energy, self.Grid_states))
            self.Grid_states[tuple(Vs)] = 1
        
        return defects_list

    
    def time_event_track(self,t,chosen_event):
        
        # Track time
        self.time.append(self.time[-1] + t)
        
        #Track events
        self.events[chosen_event] += 1
        
        
    def plot_particles(self):
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        x,y,z =  np.nonzero(self.Grid_states)
        x = x * self.steps[0]
        y = y * self.steps[1]
        z = z * self.steps[2]

        ax.scatter3D(x, y, z, c='blue', marker='o')
        
        x_lim = self.Grid_states.shape[0] * self.steps[0]
        y_lim = self.Grid_states.shape[1] * self.steps[1]
        z_lim = self.Grid_states.shape[2] * self.steps[2]
        
        ax.set_xlim(0, x_lim)
        ax.set_ylim(0, y_lim)
        ax.set_zlim(0, z_lim)
        
        ax.set_xlabel('x-axis (nm)')
        ax.set_ylabel('y-axis (nm)')
        ax.set_zlabel('z-axis (nm)')
        
        
        """
         -------------  Electrodes ---------------
        """
        left_electrode_vertices = [
        [0, self.electrodes[0] * self.steps[0], self.electrodes[0] * self.steps[0], 0, 0],
        [0, 0, y_lim, y_lim, 0],
        [z_lim, z_lim, z_lim, z_lim, z_lim]
    ]

        right_electrode_vertices = [
        [self.electrodes[1] * self.steps[0], x_lim, x_lim, self.electrodes[1] * self.steps[0], self.electrodes[1] * self.steps[0]],
        [0, 0, y_lim, y_lim, 0],
        [z_lim, z_lim, z_lim, z_lim, z_lim]
    ]
        
        electrode_vertices = [left_electrode_vertices, right_electrode_vertices]
    
        for vertices in electrode_vertices:
            # Create the 3D collection of polygons for the rectangular prism
            poly3d = Poly3DCollection([list(zip(*vertices))], facecolors='b', alpha=0.5)
            poly3d.set_edgecolor('k')  # Set edge color for better visibility
    
            # Add the collection to the plot
            ax.add_collection3d(poly3d)
        
        
        plt.show()
        
        
        
    """
    --------------------------------------------------------------------------
    Solve Poisson equation
    --------------------------------------------------------------------------
    """
    def SolvePotentialAndField(self,V):
        
        nx, ny, nz = self.Grid_states.shape

        
        num = round(nx*ny*nz) ** (1/3)
        
        u = self.u.copy()
        u1 = self.u
        
        u = self.dirichtlet(V,u) # Dirichtlet boundary conditions
        ro_charge = self.charge_distribution() # Charge distribution

        count = 0
        tol = 0.01
        
        while count < 100:
            u2 = u.copy()
            
            if count < 40:
                NN = int(60*num)
            elif (count >= 40) and (count < 70):
                NN = int(30*num)
            elif (count >= 70) and (count < 90):
                NN = int(15*num)
            elif (count >= 90) and (count < 100):
                NN = int(5*num)
                
            for _ in range (NN):

                u = self.EcPoisson(u,u1,ro_charge)
                u = self.neuman(u,u1)
                u1 = u

            count = self.check_convergence(u,u2,tol,nx,ny,nz)
            
        self.u = self.corners(u)
        
        self.Electric_field()
            
        
        
    def dirichtlet(self,V,u):
        
        u[self.electrodes[1]:,:,-1] = V # Electrode on the left grounded
        u[:,:,0] = 0 # Bottom of the active layer grounded

        u[:self.electrodes[0],:,-1] = 0 # Electrode on the right with voltage V
        

        return u
    
    def charge_distribution(self):
        
        
        qe=1.602176565e-19 # Charge of an electron (C)
        e0=8.8541878176e-12 # Vacuum permittivity (F/m)
        vol = np.prod(self.steps) * (1E-9)
        
        scale_factor = self.steps[0]*self.steps[1]
        density_particle = scale_factor* self.screening * self.q * qe / (vol*e0*self.er)
        
        
        return self.Grid_states * density_particle
    
    def EcPoisson(self,u,u1,ro_charge):
        
        c1 = 1/6
        # Vectorize implementation using indexing and slicing
        # [1:-1] index exclude the first and last element -> Elements that come 
        # from boundary conditions
        u[1:-1, 1:-1, 1:-1] = c1 * (u1[2:, 1:-1, 1:-1] + u1[:-2, 1:-1, 1:-1] +
                            u1[1:-1, 2:, 1:-1] + u1[1:-1, :-2, 1:-1] +
                            u1[1:-1, 1:-1, 2:] + u1[1:-1, 1:-1, :-2] +
                            ro_charge[1:-1, 1:-1, 1:-1])
        
        return u
        
    def neuman(self,u,u1):
        
        # Neuman conditions on the lateral faces
        u[0,1:-1,1:-1] = u1[1,1:-1,1:-1] # x = 0
        u[-1,1:-1,1:-1] = u1[-2,1:-1,1:-1] # maximum x
        u[1:-1,0,1:-1] = u1[1:-1,1,1:-1] # y = 0
        u[1:-1,-1,1:-1] = u1[1:-1,-2,1:-1] # maximum y
        
        # Neumann boundary condition in the bottom layer
        #u[1:-1,1:-1,0] = u1[1:-1,1:-1,0] # maximum y

        # Z maximum - space between electrodes
        u[self.electrodes[0]:self.electrodes[1],1:-1,-1] = u1[self.electrodes[0]:self.electrodes[1],1:-1,-2]
        
        return u
        
    def check_convergence(self,u,u2,tol,nx,ny,nz):

        total = nx * ny * nz
        count = 0
        nonzero_u = u != 0
        
        diff = np.abs(u - u2)
        diff[nonzero_u] /= u[nonzero_u]
        count += np.sum(tol > diff)
        
        return 100 * count/total
    
    def corners(self,u):
        
        u[0,0,1:-1] = 0.5 * (u[0,1,1:-1] + u[1,0,1:-1])
        u[0,-1,1:-1] = 0.5 * (u[0,-2,1:-1] - u[1,-1,1:-1])
        u[-1,0,1:-1] = 0.5 * (u[-1,1,1:-1] - u[-2,0,1:-1])
        u[-1,-1,1:-1] = 0.5 * (u[-1,-2,1:-1] - u[-2,-1,1:-1])
        
        return u

    def Electric_field(self):
        
        # Introduce the steps in m
        Ex, Ey, Ez = np.gradient(self.u,self.steps[0]*1e-9,self.steps[1]*1e-9,self.steps[2]*1e-9)
        
        self.Ex = -self.local_field(Ex)
        self.Ey = -self.local_field(Ey)
        self.Ez = -self.local_field(Ez)
        
    def local_field(self,E_field):
        
           
        """
        McPherson, J., J. Y. Kim, A. Shanware, and H. Mogul. 
        "Thermochemical description of dielectric breakdown in high dielectric constant materials." 
        Applied Physics Letters 82, no. 13 (2003): 2121-2123.
        Local electric field --> E_loc = E * (2+k)/3
        
        McPherson, J. W., and H. C. Mogul. 
        "Underlying physics of the thermochemical E model in describing low-field time-dependent dielectric breakdown in SiO 2 thin films." 
        Journal of Applied Physics 84, no. 3 (1998): 1513-1523.
        Polarization factor
        
        Padovani, A., Larcher, L., Pirrotta, O., Vandelli, L., & Bersuker, G. (2015). 
        Microscopic modeling of HfO x RRAM operations: From forming to switching. 
        IEEE Transactions on electron devices, 62(6), 1998-2006.
        
        
        # Wierzbowski, J., Klein, J., Kaniber, M., Müller, K., & Finley, J. J. 
        # Polarization control in few-layer MoS2 by electric field induced symmetry breaking.
        
        Klein, Julian, Jakob Wierzbowski, Armin Regler, Jonathan Becker, Florian Heimbach, K. Muller, Michael Kaniber, and Jonathan J. Finley. 
        "Stark effect spectroscopy of mono-and few-layer MoS2." 
        Nano letters 16, no. 3 (2016): 1554-1559.
        p = 0.08 D = 0.0166552e-10 (em)
        """
        # Dipole moment --> (enm)    %p=40;
        # p=40;
        # p=15
        
        #p = 0.0166552e-10
        #p = 1.266552e-10
    
        L=1/3
        susceptibility=self.er-1
        
        local_field=E_field*(1+L*susceptibility)
        
        return local_field
    
    def density_defects(self):
        
        nx, ny, nz = self.Grid_states.shape
        density_Vs = np.zeros((nx,ny,nz))
        square_size = self.square_size
        
        vol = (square_size ** 3) * np.prod(self.steps) 
                
        for i in range(nx):
            start_i = max(0, i - square_size // 2)
            end_i = min(nx, start_i + square_size)
            
            for j in range(ny):
                start_j = max(0, j - square_size // 2)
                end_j = min(ny, start_j + square_size)
                
                for k in range(nz):
                    start_k = max(0, k - square_size // 2)
                    end_k = min(nz, start_k + square_size)
    
                
                        
                    
                    density_Vs[i,j,k] = np.sum(self.Grid_states[start_i:end_i,start_j:end_j,start_k:end_k])/vol
    
        self.density_Vs = density_Vs
        
    def Schottky_current_2(self,V):    
        
        
        n = self.n
        nm3_to_cm3 = 1E21
        n = n * nm3_to_cm3
        
        
        region_doping = round(self.w * 1E9 / self.steps[2]) # Number of layers

        defects_at_electrodes = np.sum(self.Grid_states[:(self.electrodes[0] + region_doping + 1),:,-region_doping:])
       
        
        n1 = defects_at_electrodes/self.vol_electrode
        n1 = n1 * nm3_to_cm3
        
        delta_n = max(n1-n,0)
        
        
        # w in m --> change to cm: w*100 (cm)
        doping_term = self.e_er * np.sqrt(self.w *100 * delta_n/(4*np.pi))
        
        bias_term = np.sqrt(self.e_er/(4*np.pi)) * (2*self.e_er*n * (self.phi_b0 + self.A*np.abs(V)))**(1/4)

        phi_b = self.phi_b0 - doping_term + bias_term
        
        
        """
            IN CURRENT DENSITY (1 - np.exp(-abs(V)/self.kb_T)) SHOULD IT BE ABSOLUTE VALUE????
            Not corret --> I need to include SBH1 + SBH2 + resistance of the channel (MoS2)
            Two Schottky barriers: Metal to MoS2 + MoS2 + MoS2 to metal
            
            Coupled equations, with a voltage drop in the two contacts and the channel
            
            Self-regulated by the rectifying behavior of the SBH
            
            Zhou, Hangbo, Viacheslav Sorkin, Shuai Chen, ZhiGen Yu, Kah‐Wee Ang, and Yong‐Wei Zhang. 
            "Design‐Dependent Switching Mechanisms of Schottky‐Barrier‐Modulated Memristors based on 2D Semiconductor." 
            Advanced Electronic Materials (2023): 2201252.
        """
        current_density = self.A0 * self.T**2 * np.exp(- phi_b / (self.kb_T)) * (1 - np.exp(-abs(V)/self.kb_T))

        
        current = current_density * self.cross_sectional_area
        #print(phi_b,doping_term,bias_term,current)

        
        return current
    
    def Schottky_current(self,V):
        
        field_interface_1 = np.mean(self.Ez[:self.electrodes[0],:,-4:])
        field_interface_2 = np.mean(self.Ez[self.electrodes[1]:,:,-4:])

        
        # Due to the electric field produce by the doping and the bias
        image_force_lowering_1 = np.sqrt(self.e_er * abs(field_interface_1) / (4 * np.pi))
        image_force_lowering_2 = np.sqrt(self.e_er * abs(field_interface_2) / (4 * np.pi))

        phi_b1 = self.phi_b0 - image_force_lowering_1
        exp_phi1 = np.exp(phi_b1/self.kb_T)

        phi_b2 = self.phi_b0 - image_force_lowering_2
        exp_phi2 = np.exp(phi_b2/self.kb_T)
        
        """
        Bisection method to solve the non-linear equation
        Device structure: 
            SBH1 -> Metal-TMDC
            Channel resistance
            SBH2 -> TMDC-Metal
            
        Forward bias -> SBH2 -> TMDC-Metal dominates the total resistance
        Reverse bias -> SBH1 -> Metal-TMDC dominates the total resistance
        """
        I0 = self.cross_sectional_area * self.A0 * self.T**2 
        Vt = self.ideality_factor * self.kb_T
        
        Imin = 0

        if V>0:
            I1 = I0 * (1/exp_phi1) * (np.exp(V/Vt) - 1)
            I2 = -I0 * (1/exp_phi2) * (np.exp(-V/Vt) - 1)
            I3 = V/self.R
            Imax = min(I1, I2, I3) 
            fmax = Vt * np.log((Imax/I0) * exp_phi1 + 1) - Vt * np.log(1 - (Imax/I0) * exp_phi2) + Imax*self.R - V

        else:
            I1 = -I0 * (1/exp_phi1) * (np.exp(-V/Vt) - 1)
            I2 = I0 * (1/exp_phi2) * (np.exp(V/Vt) - 1)
            I3 = V/self.R
            Imax = max(I1, I2, I3) 
            fmax = -Vt * np.log((-Imax/I0) * exp_phi1 + 1) + Vt * np.log(1 + (Imax/I0) * exp_phi2) + Imax*self.R - V

        I = (Imax + Imin)/2
        f = 1
        
        while abs(f) >= self.tol:
            
            if V > 0:
                f = Vt * np.log((I/I0) * exp_phi1 + 1) - Vt * np.log(1 - (I/I0) * exp_phi2) + I*self.R - V
            else:
                f = -Vt * np.log(1 - (I/I0) * exp_phi1) + Vt * np.log(1 + (I/I0) * exp_phi2) + I*self.R - V
            
            if fmax > 0:
                if f > 0:
                    Imax = I
                elif f < 0:
                    Imin = I
                else:
                    break
        
            else:
                if f < 0:
                    Imax = I
                elif f > 0:
                    Imin = I
                else:
                    break
        
            I = (Imax + Imin)/2
        
        return I

    