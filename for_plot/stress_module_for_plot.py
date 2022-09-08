import numpy as np
import pandas as pd
import random
import math 
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
import csv

from mesh_distribution import MffProcessing, MeshApertureDistribution     

# ---- Stress loading module ----
class StressToAperture(MeshApertureDistribution):
    def __init__(self, node_df, ele_df, jrc_func, jcs, phi_r, Ks):
        self.ele_df = ele_df
        self.node_df = node_df
        self.jrc_func = jrc_func
        self.JCS = jcs
        self.phi_r = phi_r
        self.Ks = Ks
        
    def bondary_stress(self, stress_x, stress_y, stress_z):
        self.model_stress = [
                [(1, 0, 0), stress_x],
                [(0, 1, 0), stress_y],
                [(0, 0, 1), stress_z]
                ]
    
    def ele_normal_vector(self, ele_node):
        """ Get normal vector of element """
        tem_vector = []
        for i in range(np.shape(ele_node)[0] - 1):
            tem_vector.append(ele_node[i + 1] - ele_node[0])

        normal_vector = np.cross(tem_vector[0], tem_vector[1])
        
        return normal_vector
    
    def ele_stress(self, ele_center, normal_vector):
        """ Get stress on element """
        model_stress = self.model_stress
        
        # Get the stress of element from the boundary conditions
        # normal stress: calculate each component of stress onto the element center
        # shear stress: get the composition of model stress (vetor of composition of model stress)
        total_normal_stress = list()
        compo_vector = list()
        for stress_idx in range(len(model_stress)):
            if normal_vector[stress_idx] < 0:
                normal_vector = [normal_vector * -1 for normal_vector in normal_vector]
            # caculate cosine between normal vector and boundary stress
            cosine_demo = (
                        (np.dot(model_stress[stress_idx][0], model_stress[stress_idx][0]))**0.5 \
                        * (np.dot(normal_vector, normal_vector))**0.5
            )
            cosine_normal_stress = np.dot(model_stress[stress_idx][0], normal_vector) / cosine_demo

            if stress_idx in [0, 1]:
                model_ele_stress = (model_stress[stress_idx][1][0] * ele_center[2]) + model_stress[stress_idx][1][1]
            else:
                model_ele_stress = model_stress[stress_idx][1][1]
            total_normal_stress.append(np.array(model_ele_stress) * cosine_normal_stress)
            compo_vector.append(model_ele_stress)
           
        # normal stress: sum of the stress conduct on the element
        # shear stress: calulate the shear stress from the composition stress (parallel to the plane) 
        if normal_vector[2] < 0:
            normal_vector = [normal_vector * -1 for normal_vector in normal_vector]
            
        vector_magnitude = np.linalg.norm(compo_vector)
        tem_demo = (
            (np.dot(compo_vector, compo_vector))**0.5 \
            * (np.dot(normal_vector, normal_vector))**0.5
        )
        cosine_normal_compo = np.dot(compo_vector, normal_vector) / tem_demo
        sine_normal_compo = (1 - cosine_normal_compo)**0.5
        
        # normal stress \ shear stress
        total_normal_stress = np.sum(total_normal_stress)
        total_shear_stress = vector_magnitude * sine_normal_compo
        
        return total_normal_stress, total_shear_stress

    def cal_aperture_barton(self, ele_aper, total_normal_stress, total_shear_stress, JRC, JCS, phi_r, Ks):
        """ Normal stress: Barton & Bandis, 1983
            Shear stress: Olsson & Barton, 2001
        """ 
        # --- Part1: Normal stress
        E = ((JRC**2.5) * ele_aper * 1000000)**0.5          # unit: um
        E = E / 1000                                        # unit: mm 
        K = 0.0178 * (JCS / E) + 1.748*JRC - 7.155          # unit: Gpa ?
        a = 1 / K
        umax = -0.3 - (0.006*JRC) + 2.24*((JCS/E)**-0.25)   # unit: mm
        b = 1 / (umax * K)
        
        un = (total_normal_stress * a) / (total_normal_stress * b + 1)  # un unit: mm
        if un > umax:
            un = umax
        
        # print(ele_aper)
            
        # --- Part2: Shear stress
        # From mobilized JRC turn into mechanic aperture
        # Dilation equals to: shear displacement * tan(d_mob)
        # shear stress
        tau_peak = total_normal_stress\
                * math.tan((JRC * math.log10(JCS/total_normal_stress) + phi_r) * math.pi/180)
        tau_mob = total_shear_stress
        if tau_mob > tau_peak:
            tau_mob = tau_peak
            
        # JRC
        JRC_mob = (math.degrees(math.atan(tau_mob / total_normal_stress)) - phi_r)\
                / math.log10(JCS / total_normal_stress)
        
        # shear displacement
        shear_dis_peak = 10  # Temp unit(mm)
        shear_dis = tau_mob / Ks
        if shear_dis > shear_dis_peak:
            shear_dis = shear_dis_peak
          
        # normal displacement
        M = 2
        d_mob = 1 / M * JRC_mob * math.log10(JCS / total_normal_stress)
        nor_dis = shear_dis * math.tan((d_mob) / 180 * math.pi)
        if tau_mob == tau_peak:
            nor_dis = tau_peak / Ks * math.tan((d_mob) / 180 * math.pi)
        
        # --- Part3: calculate mechanic aperture      
        # normal stress closure
        cal_ele_m_aper = E - un
        # shear dilation
        cal_ele_m_aper += nor_dis
        
        if cal_ele_m_aper < 0:
            cal_ele_m_aper = 1e-8
        
        # ---- Part4: tranfer to hydraulic aperture (m)
        self.m_aperture = cal_ele_m_aper
        cal_hydro_aperture = (((cal_ele_m_aper * 1000)**2) / (JRC**2.5)) * (10**-6)
        
        self.un = un
        self.nor_dis = nor_dis
        self.shear = tau_mob
        self.shear_dis = shear_dis
        self.hydro_aperture = cal_hydro_aperture
        self.test = cal_ele_m_aper

        return cal_hydro_aperture
    
    def tracelength_jrc(self, fracture_index, trace_length):
        """ determine the JRC correspond to trace length """
        mean_func = self.jrc_func[fracture_index-1]
        jrc_mean = (mean_func[0]**2) * trace_length + mean_func[1]*trace_length + mean_func[2]
        
        return jrc_mean
    
    def single_ele_stress_apply(self, fracture_index, ele_idx, jrc_plane):
        jrc = jrc_plane
        jcs = self.JCS[fracture_index - 1]
        phi_r = self.phi_r[fracture_index - 1]
        Ks = self.Ks[fracture_index - 1]
        
        ele_node = np.empty([3,3])
        ele_aper = self.ele_df['Apert'][ele_idx]

        # Store coordinate of each node in element
        for node_idx, node in enumerate(self.ele_df.loc[ele_idx][0 : 3]):
            for cor_idx, cor in enumerate(self.node_df.loc[node][0 : 3]):
                ele_node[node_idx, cor_idx] = cor
        
        # Store center of element 
        ele_center = []
        for i in range(np.shape(ele_node)[1]):
            tem_center = 0
            for j in range(np.shape(ele_node)[0]):
                tem_center += ele_node[j, i]
            ele_center.append(1/3 * (tem_center))

        # Noramal vector
        normal_vector =  self.ele_normal_vector(ele_node)

        # Stress of element
        total_normal_stress, total_shear_stress = self.ele_stress(ele_center, normal_vector)
        
        # Aperture, tansmisity change of element
        ele_hydro_aper = self.cal_aperture_barton(
            ele_aper, 
            total_normal_stress, 
            total_shear_stress,
            jrc, 
            jcs,
            phi_r,
            Ks
            )

        return ele_hydro_aper
    
    def main_process(self):
        """
        1. Read though the df of element
            1.1. Find the mesh belong to one single plane
            1.2. Merge into one polygon
            1.3. Determine the trace length and JRC of plane
            1.4. Apply stress with corresponding JRC
        """
        # Read though the df of element
        fracture_index = 1 
        fracture_set_max = max(self.ele_df['Set#'])
        output_df = self.ele_df.copy()
        while True:
            # 1.1. Find the mesh belong to one single plane
            tem_node_list = self.collect_node_data(fracture_index)
            
            # 1.2. Merge into one polygon
            merge_polygon = self.merge_meshes(tem_node_list)
        
            # 1.3. Determine the trace length and JRC of plane
            trace_length = self.polygon_tracelength(merge_polygon)
            jrc_plane = self.tracelength_jrc(fracture_index, trace_length)
            
            # 1.4. Apply stress with corresponding JRC
            # And the corresponding trans, stor
            for ele_idx in self.ele_df[self.ele_df['Set#'] == fracture_index].index:
                tem_hyd_aperture = self.single_ele_stress_apply(fracture_index, ele_idx, jrc_plane)
                
                output_df['Apert'][ele_idx] = tem_hyd_aperture
                output_df['Trans'][ele_idx] = 817500 * tem_hyd_aperture**3
            fracture_index += 1
            if fracture_index > fracture_set_max:
                return output_df
    
if __name__ == '__main__':
    #path = '/Users/weitim/Desktop/研究/研究主軸/FracMan_module/stress_on_aperture_distriution/'
    file_name = 'FlowDefinition_only_H.mff'
    node_df, ele_df = MffProcessing.path_input(file_name)
    
    # fracture strength parameters
    # trace_jrc = np.array([
    #     [0.463, -1.1643, 12.456],
    #     [0.463, -1.1643, 12.456],
    # ])
    trace_jrc = np.array([
        [0, 0, 10],
        [0, 0, 10],
    ])
    jcs = [100, 100]
    phi_r = [30, 30]
    Ks = [2, 2]
    
    test_model_stress = StressToAperture(node_df, ele_df, trace_jrc, jcs, phi_r, Ks)
    
    
    # stress condition
    z_list = np.linspace(0, 10, 20)
    
    # stress_x = [0, 10]
    # stress_y = [0, 10]
    # stress_z = [0, 10]
    # boundary = [stress_x, stress_y, stress_z]
    # test_model_stress.bondary_stress(boundary[0], boundary[1], boundary[2])
    # test_model_stress.main_process()
    nor_dis_list = []
    shear_dis_list = []
    shear_list = []
    ap_list = []
    test_list = []
    for z in z_list: 
        stress_x = [0, z]
        stress_y = [0, z]
        stress_z = [0, 5]
        boundary = [stress_x, stress_y, stress_z]
        test_model_stress.bondary_stress(boundary[0], boundary[1], boundary[2])
        test_model_stress.main_process()
        shear_list.append(test_model_stress.shear)
        shear_dis_list.append(test_model_stress.shear_dis)
        nor_dis_list.append(test_model_stress.nor_dis)
        ap_list.append(test_model_stress.hydro_aperture)
        test_list.append(test_model_stress.test)

    
    fig, ax = plt.subplots(2,1)
    ax[0].plot(shear_dis_list, shear_list, 'b-.o')  
    ax[0].set_ylabel('Shear stress (MPa)')
    ax[0].set_xlabel('Shear displacement (mm)')  
    ax[0].set_title('Shear stress v.s. Shear displacement')
    
    # fig, ax = plt.subplots()
    ax[1].plot(shear_dis_list, nor_dis_list, 'b-.o')  
    ax[1].set_ylabel('Dilation (mm)')
    ax[1].set_xlabel('Shear displacement (mm)')  
    ax[1].set_title('Dilation v.s. Shear displacement')
    
    # ax[2].plot(shear_list, test_list, 'b-.o')
    # ax[2].set_ylabel('Mechanical aperture (mm)')
    # ax[2].set_xlabel('Shear stress (Mpa)')  
    # ax[2].set_title('Mechanical aperture v.s. Shear stress')
    
    # ax[3].plot(shear_list, test_list, 'b-.^')
    # ax[3].set_ylabel('Hydraulic aperture (mm)')
    # ax[3].set_xlabel('Shear stress (Mpa)')  
    # ax[3].set_title(r"$\sigma_{n}$")
    # ax.set_ylim([-0.1, 0.2])
    fig.tight_layout()
    plt.show() 
    # print(ap_list)
    
    # fig, ax = plt.subplots()
    # ax.plot(z_list, ap_list, 'b-.o')  
    # ax.set_xlabel('Normal stress (Mpa)')
    # ax.set_ylabel('Mechanic aperture (mm)')  
    # ax.set_title('Mechanic aperture v.s. Normal stress')
    # plt.show()    
    