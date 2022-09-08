from cgi import test
from turtle import color
import numpy as np
import pandas as pd
import random
import math 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
import csv

from mesh_distribution import MffProcessing, MeshApertureDistribution     

# ---- Stress loading module ----
class StressToAperture(MeshApertureDistribution):
    def __init__(self, node_df, ele_df, jrc_func, jcs, phi_r):
        self.ele_df = ele_df
        self.node_df = node_df
        self.jrc_func = jrc_func
        self.JCS = jcs
        self.phi_r = phi_r
        
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
        normal_vector /= np.linalg.norm(normal_vector)
        
        return normal_vector
    
    def ele_stress(self, ele_centroid, normal_vector):
        """ Get stress on element """
        
        model_stress = self.model_stress
        # Part1: Use Cauchy's formula (F = \sigma \cdot N) to calculate the traction from boundary stress
        # Cauchy stress tensor (variation: sigma)
        sigma = np.zeros([3, 3])
        for stress_idx in range(len(model_stress)):
            # Get the element boundary stress from boundary condistion
            # stress relate with z
            ele_boundary_stress = (model_stress[stress_idx][1][0] * ele_centroid[2]) + model_stress[stress_idx][1][1]
            sigma[stress_idx, stress_idx] = ele_boundary_stress

        # Part2: Calculate the traction(variation: traction) of fracture element from element Cauchy stress tensor
        traction = np.dot(sigma, normal_vector)
        
        # Part3: Calculate the normal stress & shear stress from traction of fracture element 
        element_normal_stress = traction.dot(normal_vector)
        traction_normal_vector_cosine = traction.dot(normal_vector) \
                                        / (np.linalg.norm(traction) * np.linalg.norm(normal_vector))
        if traction_normal_vector_cosine > 1:
            traction_normal_vector_cosine = 1
            
        traction_normal_vector_theta = np.arccos(traction_normal_vector_cosine)
        traction_normal_vector_sine = np.sin(traction_normal_vector_theta)
        element_shear_stress = np.linalg.norm(traction) * traction_normal_vector_sine
        
        return element_normal_stress, element_shear_stress
       
    def cal_aperture_barton(self, ele_aper, ele_normal_stress, ele_shear_stress, JRC, JCS, phi_r):
        """ Normal stress: Barton & Bandis, 1983
            Shear stress: Olsson & Barton, 2001
        """ 

        E = ((JRC**2.5) * ele_aper * 1e6)**0.5          # unit: um
        E = E / 1000                                        # unit: mm 
        # Part1: Normal stress
        un = self.normal_stress_barton(E, ele_normal_stress, JRC, JCS)

        # Part2: shear stress
        nor_dis = self.shear_stress_barton(ele_normal_stress, ele_shear_stress, JRC, JCS, phi_r)
  
        # Part3: calculate mechanic aperture      
        #   normal stress closure [mm]
        cal_ele_m_aper = E - un
        #   shear dilation [mm]
        cal_ele_m_aper += nor_dis
        
        if cal_ele_m_aper < 0:
            cal_ele_m_aper = 6e-4
        
        # Part4: tranfer to hydraulic aperture (m)
        cal_hydro_aperture = ((cal_ele_m_aper * 1000)**2) / (JRC**2.5)  # unit: um
        cal_hydro_aperture /= (10**6)   # unit: m

        return cal_hydro_aperture

    def normal_stress_barton(self, E, ele_normal_stress, JRC, JCS):
        """ Express fracture closure behavior, following Barton et al., 1983 
        Args:
            E (float): Mechanical aperture of element [mm]
            ele_normal_stress (float): normal stress conduct on element
            JRC (float): JRC of fracture
            JCS (float): JCS of fracture
        Returns:
            float(un): closure of fracture under normal stress [mm]
        """
        K = 0.0178 * (JCS / E) + 1.748*JRC - 7.155          # unit: Gpa ?
        a = 1 / K
        umax = -0.3 - (0.006*JRC) + 2.24*((JCS/E)**-0.25)   # unit: mm
        b = 1 / (umax * K)
        un = (ele_normal_stress * a) / (ele_normal_stress * b + 1)  # un unit: mm
        if un > umax:
            return umax
        else:
            return un
            
    def shear_stress_barton(self, ele_normal_stress, ele_shear_stress, JRC, JCS, phi_r):
        """ Express fracture dilation behavior under shear stress, following Barton et al., 1982 
        Args:
            ele_normal_stress (float): _description_
            ele_shear_stress (float): _description_
            JRC (float): _description_
            JCS (float): _description_
            phi_r (float): _description_

        Returns:
            float: _description_
        """
        if JCS == ele_normal_stress:
            ele_normal_stress -= 0.001
        i_theta = JRC * math.log10(JCS / ele_normal_stress)
        jrc_table = (-(phi_r / i_theta), 0, 0.75, 1, 0.85, 0.7, 0.5, 0)
        dis_table = (0, 0.3, 0.6, 1, 2, 4, 10, 100)  
        
        tau_peak = ele_normal_stress\
                * math.tan((JRC * math.log10(JCS / ele_normal_stress) + phi_r) * math.pi/180)
        tau_mob = ele_shear_stress
        if ele_shear_stress > tau_peak:
            tau_mob = tau_peak
            
        # JRC
        JRC_mob = (math.degrees(math.atan(tau_mob / ele_normal_stress)) - phi_r)\
                / math.log10(JCS / ele_normal_stress)
        JRC_peak = JRC
        if JRC_mob > JRC_peak:
            JRC_mob = JRC_peak
        jrc_mob_div_jrc_peak = JRC_mob/ JRC_peak

        # shear displacement
        shear_dis_peak = 1000/500 * ((JCS/1000) ** 0.33)  # Temp unit(mm)
        # shear_dis_peak = 2
        
        # Use the table from Barton et al., 1982 to calculate the shear displacement
        for i in range(len(jrc_table)-1):
            if jrc_mob_div_jrc_peak < 1:
                if (jrc_mob_div_jrc_peak >= jrc_table[i]) and (jrc_mob_div_jrc_peak < jrc_table[i+1]):
                    jrc_ratio = (jrc_mob_div_jrc_peak - jrc_table[i]) / (jrc_table[i+1] - jrc_table[i])
                    dis_ratio = jrc_ratio * (dis_table[i+1] - dis_table[i]) + dis_table[i]
                    shear_dis = dis_ratio * shear_dis_peak
                elif jrc_mob_div_jrc_peak < jrc_table[0]:
                    shear_dis = 0        
            elif jrc_mob_div_jrc_peak >= 1:
                if (jrc_mob_div_jrc_peak <= jrc_table[i]) and (jrc_mob_div_jrc_peak > jrc_table[i+1]):
                    jrc_ratio = (jrc_table[i] - jrc_mob_div_jrc_peak) / (jrc_table[i] - jrc_table[i+1])
                    dis_ratio = jrc_ratio * (dis_table[i] - dis_table[i+1]) + dis_table[i]
                    shear_dis = dis_ratio * shear_dis_peak
                elif jrc_mob_div_jrc_peak <= jrc_table[-1]:
                    shear_dis = 100 * shear_dis_peak
                
        # normal displacement
        M = 2
        d_mob = 1 / M * JRC_mob * math.log10(JCS / ele_normal_stress)
        nor_dis = shear_dis * math.tan((d_mob) / 180 * math.pi)
        
        return nor_dis
    
    def tracelength_jrc(self, fracture_index, trace_length):
        """ determine the JRC correspond to trace length """
        mean_func = self.jrc_func[fracture_index-1]
        jrc_mean = mean_func[0] * (trace_length**2) + mean_func[1] * trace_length + mean_func[2]
        return jrc_mean
    
    def element_centroid(self, ele_node):
        ''' Get the element centroid '''
        ele_center = []
        for i in range(np.shape(ele_node)[1]):
            tem_center = 0
            for j in range(np.shape(ele_node)[0]):
                tem_center += ele_node[j, i]
            ele_center.append(1/3 * (tem_center))
            
        return ele_center
        
    def single_ele_stress_apply(self, fracture_index, ele_idx, jrc_plane):
        jrc = jrc_plane
        jcs = self.JCS[fracture_index - 1]
        phi_r = self.phi_r[fracture_index - 1]
        
        ele_node = np.empty([3,3])
        ele_aper = self.ele_df['Apert'][ele_idx]

        # Store coordinate of each node in element
        for node_idx, node in enumerate(self.ele_df.loc[ele_idx][0 : 3]):
            for cor_idx, cor in enumerate(self.node_df.loc[node][0 : 3]):
                ele_node[node_idx, cor_idx] = cor
        
        # Get the centroid of element 
        ele_centroid = self.element_centroid(ele_node)

        # Noramal vector
        normal_vector =  self.ele_normal_vector(ele_node)

        # Stress of element
        element_normal_stress, element_shear_stress = self.ele_stress(ele_centroid, normal_vector)
        
        # Aperture, tansmisity change of element
        ele_hydro_aper = self.cal_aperture_barton(
            ele_aper, 
            element_normal_stress, 
            element_shear_stress,
            jrc, 
            jcs,
            phi_r
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
                output_df['Trans'][ele_idx] = 817500 * (tem_hyd_aperture**3)
            fracture_index += 1
            if fracture_index > fracture_set_max:
                return output_df
    

if __name__ == '__main__':
    pass
    # output_df = test_model_stress.main_process()
    # print(output_df)
        