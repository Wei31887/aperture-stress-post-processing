import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mesh_distribution import AperturePostProcessing, HeadBoundaryAssign, MeshApertureDistribution    
from mff_processing import MffProcessing

if __name__ == '__main__':
    # import data file
    
    path = 'D:\\wei_data\Research\\REV_determination\\REV_determination\\REV_without_stress_without_dis\\8\\8_x\\'
    #path = 'D:\\wei_data\Research\\REV_determination\\REV_determination\\REV_without_stress_without_dis\\20\\1\\20_1_x\\'
    file_name = 'boundary_8_x.mff'
    
    # import mesh data from mff.file
    node_df, ele_df = MffProcessing.path_input(path + file_name)

    # aperture distribution
    trace_aperture = np.array([
        [[-5.70E-07, 2.25E-05, 9.78E-06], [-4.46E-07, 1.40E-05, 7.20E-07]],
        [[-4.49E-07, 1.77E-05, 1.86E-05], [-1.56E-07, 5.84E-06, 2.32E-06]],
        [[-7.67E-07, 3.06E-05, 4.58E-07], [-4.33E-07, 1.16E-05, -4.72E-06]]
    ])

    aperture_dis = MeshApertureDistribution(node_df, ele_df, trace_aperture)
    cal_ele_df = aperture_dis.main_distribution()
    
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
    MffProcessing.ele_output_mff(path + file_name, 'ele.mff', cal_ele_df)
    
    
    # node 
    MffProcessing.node_output_mff('ele.mff', 'output_x_8.mff', cal_node_df_x)
    MffProcessing.node_output_mff('ele.mff', 'output_y_8.mff', cal_node_df_y)
    MffProcessing.node_output_mff('ele.mff', 'output_z_8.mff', cal_node_df_z)