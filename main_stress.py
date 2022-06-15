import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mesh_distribution import HeadBoundaryAssign
from mff_processing import MffProcessing
from stress_module import StressToAperture

if __name__ == '__main__':
    # import data file

    # path = '/Users/weitim/Desktop/研究/研究主軸/FracMan_module/stress_on_aperture_distriution/'
    path = 'D:\\wei_data\\Research\\REV_determination\\REV_determination\\REV_with_stress_with_dis\\for_stress_origin_5\\origin_mff\\'
    file_name = 'boundary_16_x.mff'

    #node, element df
    node_df, ele_df = MffProcessing.path_input(path + file_name)
    
    # ---------boundary stress
    trace_jrc = np.array([
        [0.0233, -0.7578, 10.972],
        [0.0274, -0.9533, 13.109],
        [0.0305, -1.1256, 16.561]
    ])
    jcs = [16, 11.5, 10.7]
    phi_r = [32.7, 38, 35.8]
    
    model_stress = StressToAperture(node_df, ele_df, trace_jrc, jcs, phi_r)
    stress_x = [0, 1.6]
    stress_y = [0, 1.6]
    stress_z = [0, 0.4]
    boundary = [stress_x, stress_y, stress_z]
    model_stress.bondary_stress(boundary[0], boundary[1], boundary[2])
    cal_ele_df = model_stress.main_process()
    
    # output calculated element file
    MffProcessing.ele_output_mff(path + file_name, 'ele.mff', cal_ele_df)
    
    # ----------

    # three different boundaries
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

    # node 
    MffProcessing.node_output_mff('ele.mff', 'stress_k_4_x.mff', cal_node_df_x)
    MffProcessing.node_output_mff('ele.mff', 'stress_k_4_y.mff', cal_node_df_y)
    MffProcessing.node_output_mff('ele.mff', 'stress_k_4_z.mff', cal_node_df_z)
    