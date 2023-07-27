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
from scipy.constants import elementary_charge as qe, epsilon_0 as e0




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
        self.e_er = qe/(self.er*e0)
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
        
        self.phi_b = [[],[]]
        
        
        """
        ------------------- Depletion region --------------------------------
        """
        # Positive pointing to the electrodes (z direction)
        # Make more difficult for the electron to cross the barrier
        ND = 1E16
        phi_bi = 0.2
        max_field = 100 * np.sqrt(2*self.e_er*100 * ND * phi_bi) #(V/m)
        WD = np.sqrt(2 * (phi_bi - self.kb_T)/(100 * self.e_er * ND)) *1e7
        
        x = np.linspace(self.steps[2],self.device_size[2],self.Grid_states.shape[2])
        
        self.depletion_region_field = max_field - x * max_field /WD
        """
        --------------------------------------------------------------------
        """


        
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
            Vs_ijk = rng.integers(1, np.array(grid_points)-2, size=3)

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
        #single_particle = [np.array([1,15,29]),np.array([1,14,29]),np.array([1,13,29])]
        single_particle = [np.array([1,15,29]),np.array([1,14,29]),np.array([1,13,29]),np.array([1,16,29]),np.array([1,17,29])]
        
        test_defects = [x_axis_paralallel,y_axis_paralallel,z_axis_paralallel,single_particle]
        
        for Vs in test_defects[3]:
            defects_list.append(Defects(Vs, self.Act_energy, self.Grid_states,self.q))
            self.Grid_states[tuple(Vs)] = 1
        
        return defects_list

    
    def time_event_track(self,t,chosen_event):
        
        # Track time
        self.time.append(self.time[-1] + t)
        
        #Track events
        self.events[chosen_event] += 1
        
    """
    --------------------- Plot the graphs ---------------------
    """
        
    def plot_particles(self,V,current,path,i,t):
                
        nr = 4
        nc = 4
        
        #ax = plt.subplot2grid(shape=(nr, nc), loc=(0, 0), rowspan=1, colspan=1)
        ax = plt.subplot2grid(shape=(nr, nc), loc=(0, 0), rowspan=2, colspan=3, projection='3d')
        ax.title.set_text(f'Voltage (V): {V.voltage[-1]:.2f}, cycle = {V.cycles}')
        #ax = fig.add_subplot(111, projection='3d')

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
        
        #ax.set_xlabel('x-axis (nm)')
        #ax.set_ylabel('y-axis (nm)')
        #ax.set_zlabel('z-axis (nm)')
        
        
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
            
        # I-V curve plot
        ax1 = plt.subplot2grid(shape=(nr, nc), loc=(0, 3), rowspan=1, colspan=1)
        ax1.plot(V.voltage,abs(np.array(current)))
        ax1.set_ylim(1e-11,1e-3)
        ax1.set_xlabel('Voltage (V)', color ='black')
        ax1.set_ylabel('Current (A)', color='black')
        #ax1.set_ylabel('Current (μA)', color='black')

        ax1.set_yscale('log')
        
        
        ax2 = plt.subplot2grid(shape=(nr, nc), loc=(2, 0), rowspan=2, colspan=3, projection='3d')
        stride = 6 # Adjust the stride value as needed
        scale_factor = 1.5  # Adjust the scale factor as needed
        m_to_nm = 1e-9
        x,y,z = np.mgrid[self.steps[0]:self.device_size[0]:self.steps[0] * stride,
                            self.steps[1]:self.device_size[1]:self.steps[1] * stride,
                            self.steps[2]:self.device_size[2]:self.steps[2] * stride]
        
        ax2.quiver3D(x, y, z, 
                    self.Ex[::stride,::stride,::stride] * m_to_nm, 
                    self.Ey[::stride,::stride,::stride] * m_to_nm, 
                    self.Ez[::stride,::stride,::stride] * m_to_nm, 
                    length=scale_factor)
        
        
        ax2.set_xlabel('x axis (nm)', color ='black')
        ax2.set_ylabel('y axis (nm)', color='black')
        ax2.set_zlabel('z axis (nm)', color='black')
        
        
        ax3 = plt.subplot2grid(shape=(nr, nc), loc=(2, 3), rowspan=2, colspan=1)
        ax3.plot(V.voltage,self.phi_b[0],color = 'blue',label = 'SBH_1')
        ax3.plot(V.voltage,self.phi_b[1],color = 'red',label = 'SBH_2')
        
        ax3.set_xlabel('Voltage (V)', color ='black')
        ax3.set_ylabel('SBH (eV)', color='black')
        
        ax3.legend(fontsize="5")





        if path == '':
            plt.show()
        else:
            plt.savefig(path+str(i)+'_t(s) = '+str(round(t,5))+' .png', dpi = 300)
            plt.clf()
        
        
        
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
        
        
        vol = np.prod(self.steps) * (1E-9)
        ro_charge = np.zeros(self.Grid_states.shape)
        
        scale_factor = self.steps[0]*self.steps[1]
        density_particle_screening = scale_factor* self.screening[0] * self.q * qe / (vol*e0*self.er)
        density_particle = scale_factor* self.screening[1] * self.q * qe / (vol*e0*self.er)
        
        ro_charge = self.Grid_states * density_particle
        ro_charge[self.electrodes[0]:self.electrodes[1],:,:] = self.Grid_states[self.electrodes[0]:self.electrodes[1],:,:] * density_particle_screening
        return ro_charge
    
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
        
        # There is a bug here --> Sometimes it doesn't converge
        # diff seems okay, probably the problem raise from /u[non_zero_u]
        diff[nonzero_u] /= u[nonzero_u] 

        count += np.sum(tol > diff)
        
        return 100 * count/total
    
    def corners(self,u):
        
        u[0,0,1:-1] = 0.5 * (u[0,1,1:-1] + u[1,0,1:-1])
        u[0,-1,1:-1] = 0.5 * (u[0,-2,1:-1] + u[1,-1,1:-1])
        u[-1,0,1:-1] = 0.5 * (u[-1,1,1:-1] + u[-2,0,1:-1])
        u[-1,-1,1:-1] = 0.5 * (u[-1,-2,1:-1] + u[-2,-1,1:-1])
        
        u[self.electrodes[0]:self.electrodes[1],0,-1] = 0.5 * (u[self.electrodes[0]:self.electrodes[1],0,-2] + u[self.electrodes[0]:self.electrodes[1],1,-1])
        u[self.electrodes[0]:self.electrodes[1],-1,-1] = 0.5 * (u[self.electrodes[0]:self.electrodes[1],-1,-2] + u[self.electrodes[0]:self.electrodes[1],-2,-1])

        
        return u

    def Electric_field(self):
        
        # Introduce the steps in m
        Ex, Ey, Ez = np.gradient(self.u,self.steps[0]*1e-9,self.steps[1]*1e-9,self.steps[2]*1e-9)
        
         
        self.Ex = -self.local_field(Ex)
        self.Ey = -self.local_field(Ey)
        self.Ez = -self.local_field(Ez)
        
        self.Ez[:self.electrodes[0],:,:] = self.Ez[:self.electrodes[0],:,:] + self.depletion_region_field[::-1]
        self.Ez[self.electrodes[1]:,:,:] = self.Ez[self.electrodes[1]:,:,:] + self.depletion_region_field[::-1]
        
        
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
        1 D ≈ 0.02081943 e·nm 
        
        Yang, Longlong, Xin Xie, Jingnan Yang, Mengfei Xue, Shiyao Wu, Shan Xiao, Feilong Song et al. 
        "Strong light–matter interactions between gap plasmons and two-dimensional excitons under ambient conditions in a deterministic way." 
        Nano Letters 22, no. 6 (2022): 2177-2186.
        
        p = 7.5 D
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
        doping_term = self.e_er * 100 * np.sqrt(self.w *100 * delta_n/(4*np.pi))
        
        bias_term = np.sqrt(self.e_er * 100/(4*np.pi)) * (2*self.e_er* 100 * n * (self.phi_b0 + self.A*np.abs(V)))**(1/4)

        print(doping_term,bias_term)

        phi_b = self.phi_b0 - doping_term + bias_term
        
        
        """
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
        
        method_calculate_phi = ['potential','average_field','average_phi']
        calculate_phi_used = method_calculate_phi[2]
        methods = ['Newton_raphson','Bisection']
        method_used = methods[0]
        
        
        if calculate_phi_used == 'potential':
            field_interface_1 = np.mean(self.Ez[:self.electrodes[0],:,:],(0,1)) # Average the plane xy
            field_interface_2 = np.mean(self.Ez[self.electrodes[1]:,:,:],(0,1)) # Average the plane xy
            
            x = np.linspace(self.steps[2],self.device_size[2],self.Grid_states.shape[2])
    
            PE_1 = np.argmax(- self.e_er * qe / (16*np.pi * x * 1e-9) - qe * x * 1e-9 * abs(field_interface_1[::-1]))
            PE_2 = np.argmax(- self.e_er * qe / (16*np.pi * x * 1e-9) - qe * x * 1e-9 * abs(field_interface_2[::-1]))
            
            # Due to the electric field produce by the doping and the bias
            image_force_lowering_1 = np.sqrt(self.e_er * abs(field_interface_1[PE_1]) / (4 * np.pi))
            image_force_lowering_2 = np.sqrt(self.e_er * abs(field_interface_2[PE_2]) / (4 * np.pi))
            
            phi_b1 = self.phi_b0 - image_force_lowering_1
            phi_b2 = self.phi_b0 - image_force_lowering_2
        
        elif calculate_phi_used == 'average_field':
            field_interface_1 = np.mean(self.Ez[:self.electrodes[0],:,-6:],(0,1)) # Average the plane xy
            field_interface_2 = np.mean(self.Ez[self.electrodes[1]:,:,-6:],(0,1)) # Average the plane xy
            
            
            image_force_lowering_1 = np.sqrt(self.e_er * max(abs(field_interface_1)) / (4 * np.pi))
            image_force_lowering_2 = np.sqrt(self.e_er * max(abs(field_interface_2)) / (4 * np.pi))
            
            phi_b1 = self.phi_b0 - image_force_lowering_1
            phi_b2 = self.phi_b0 - image_force_lowering_2
            
        elif calculate_phi_used == 'average_phi':
            
            image_force_lowering_1 = np.sqrt(self.e_er * abs(self.Ez[:self.electrodes[0],:,-1]) / (4 * np.pi))
            image_force_lowering_2 = np.sqrt(self.e_er * abs(self.Ez[self.electrodes[1]:,:,-1]) / (4 * np.pi))
            #print(max(abs(field_interface_1))*1e-9,max(abs(field_interface_2))*1e-9)
            
            phi_b = np.zeros(self.Grid_states.shape[:2]) + self.phi_b0
            
            phi_b[:self.electrodes[0],:] -= image_force_lowering_1
            phi_b[self.electrodes[1]:,:] -= image_force_lowering_2 
            
            phi_b1 = np.mean(phi_b[:self.electrodes[0],:])
            phi_b2 = np.mean(phi_b[self.electrodes[1]:,:])

        

            
        exp_phi1 = np.exp(phi_b1/self.kb_T)
        exp_phi2 = np.exp(phi_b2/self.kb_T)
        
        
        self.phi_b[0].append(phi_b1)
        self.phi_b[1].append(phi_b2)
        
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
        f = 1


        if method_used == methods[0]:
            
            
            if V >= 0:
                I1 = -I0 * (1/exp_phi1) * (np.exp(-V/Vt) - 1)
                I2 = I0 * (1/exp_phi2) * (np.exp(V/Vt) - 1)
                I3 = V/self.R
                I = min(I1,I2,I3)
                
                if abs(I) < 1e-15:
                    return 0
                elif (I/I0) * exp_phi1 >= 1:
                    return I1

            elif V < 0:
                I1 = I0 * (1/exp_phi1) * (np.exp(abs(V)/Vt) - 1)
                I2 = -I0 * (1/exp_phi2) * (np.exp(-abs(V)/Vt) - 1)
                I3 = abs(V)/self.R
                I = min(I1,I2,I3)
                
                if abs(I) < 1e-15:
                    return 0
                elif (I/I0) * exp_phi2 >= 1:
                    return I2
                
                
            while (abs(f) >= self.tol * abs(V)):
                    
                if V > 0:
                    

                    if (I/I0) * exp_phi1 >= 0.99999999:
                        return I 
                    
                    f = Vt * np.log((I/I0) * exp_phi2 + 1) - Vt * np.log(1 - (I/I0) * exp_phi1) + I*self.R - V
                    df = (Vt * exp_phi2) / (I * exp_phi2 + I0) + Vt * exp_phi1 / (I0 - I * exp_phi1) + self.R
                    
                else:
                    
                    if (I/I0) * exp_phi2 >= 0.99999999:
                        return I 
                    
                    f = -Vt * np.log(1 - (I/I0) * exp_phi2) + Vt * np.log(1 + (I/I0) * exp_phi1) + I*self.R - abs(V)
                    df = (Vt * exp_phi2) / (I0 - I * exp_phi2) + Vt * exp_phi1 / (I0 + I * exp_phi1) + self.R


                I=I-f/df


        elif method_used == methods[1]:
 
    
            
            Vt = self.ideality_factor * self.kb_T
    
            Imin = 0
    
            if V>0:
                I1 = -I0 * (1/exp_phi1) * (np.exp(-V/Vt) - 1)
                I2 = I0 * (1/exp_phi2) * (np.exp(V/Vt) - 1)
                I3 = V/self.R
                Imax = min(I1, I2, I3) 
                #fmax = Vt * np.log((Imax/I0) * exp_phi1 + 1) - Vt * np.log(1 - (Imax/I0) * exp_phi2) + Imax*self.R - V
                fmax = Vt * np.log((Imax/I0) * exp_phi2 + 1) - Vt * np.log(1 - (Imax/I0) * exp_phi1) + Imax*self.R - V
    
            else:
                I1 = I0 * (1/exp_phi1) * (np.exp(V/Vt) - 1)
                I2 = -I0 * (1/exp_phi2) * (np.exp(-V/Vt) - 1)
                I3 = V/self.R
                Imax = max(I1, I2, I3) 
                #fmax = -Vt * np.log((-Imax/I0) * exp_phi1 + 1) + Vt * np.log(1 + (Imax/I0) * exp_phi2) + Imax*self.R - V
                fmax = -Vt * np.log((-Imax/I0) * exp_phi2 + 1) + Vt * np.log(1 + (Imax/I0) * exp_phi1) + Imax*self.R - V
    
            I = (Imax + Imin)/2
            
            while abs(f) >= self.tol:
                
                if V > 0:
                    #f = Vt * np.log((I/I0) * exp_phi1 + 1) - Vt * np.log(1 - (I/I0) * exp_phi2) + I*self.R - V
                    f = Vt * np.log((I/I0) * exp_phi2 + 1) - Vt * np.log(1 - (I/I0) * exp_phi1) + I*self.R - V
                else:
                    #f = -Vt * np.log(1 - (I/I0) * exp_phi1) + Vt * np.log(1 + (I/I0) * exp_phi2) + I*self.R - V
                    f = -Vt * np.log(1 - (I/I0) * exp_phi2) + Vt * np.log(1 + (I/I0) * exp_phi1) + I*self.R - V
                
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
        
        if V > 0:
            return I
        else:
            return -I

    