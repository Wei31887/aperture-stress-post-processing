import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mesh_distribution import AperturePostProcessing, HeadBoundaryAssign, MeshApertureDistribution    
from mff_processing import MffProcessing

if __name__ == '__main__':
    # import data file
    
    #path = 'D:\\wei_data\Research\\REV_determination\\REV_determination\\REV_without_stress_without_dis\\8\\8_x\\'
    path = 'D:\\wei_data\Research\\REV_determination\\REV_determination\\REV_without_stress_without_dis\\16\\5\\16_5_x\\'
    file_name = 'boundary_16_x.mff'
    
    # import mesh data from mff.file
    node_df, ele_df = MffProcessing.path_input(path + file_name)
    
    # head boundary
    region = np.array([10, 10, 10, 20, 20, 20])
    #   x
    constant_10 = [0, 0, 0, 20]
    constant_0 = [0, 0, 0, 0]
    var_x_1 = [1, 0, 0, 0]
    boundary = np.array([
        constant_10,
        constant_0,
        var_x_1
    ])
    boundary_type = 'x'

    boundary_x_assign = HeadBoundaryAssign(node_df, boundary, boundary_type, region)
    cal_node_df_x = boundary_x_assign.main_assign(err=1e-3)

    #   y
    constant_10 = [0, 0, 0, 20]
    constant_0 = [0, 0, 0, 0]
    var_y_1 = [0, 1, 0, 0]
    boundary = np.array([
        constant_10,
        constant_0,
        var_y_1
    ])
    boundary_type = 'y'

    boundary_y_assign = HeadBoundaryAssign(node_df, boundary, boundary_type, region)
    cal_node_df_y = boundary_y_assign.main_assign(err=1e-3)

    #   z
    constant_10 = [0, 0, 0, 20]
    constant_0 = [0, 0, 0, 0]
    var_z_1 = [0, 0, 1, 0]
    boundary = np.array([
        constant_10,
        constant_0,
        var_z_1
    ])
    boundary_type = 'z'

    boundary_z_assign = HeadBoundaryAssign(node_df, boundary, boundary_type, region)
    cal_node_df_z = boundary_z_assign.main_assign(err=1e-3)

    # output mff.file of three boundaries
    # element
    MffProcessing.ele_output_mff(path + file_name, 'no_dis_ele.mff', ele_df)
    
    
    # node 
    MffProcessing.node_output_mff('no_dis_ele.mff', 'nodis_x_16_5.mff', cal_node_df_x)
    MffProcessing.node_output_mff('no_dis_ele.mff', 'nodis_y_16_5.mff', cal_node_df_y)
    MffProcessing.node_output_mff('no_dis_ele.mff', 'nodis_z_16_5.mff', cal_node_df_z)