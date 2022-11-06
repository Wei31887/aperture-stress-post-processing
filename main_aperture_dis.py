import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mesh_distribution import AperturePostProcessing, HeadBoundaryAssign, MeshApertureDistribution    
from mff_processing import MffProcessing

if __name__ == '__main__':
    # 輸入路徑
    
    path = "/Users/weitim/Desktop/Research/organized_modules/aperture-postprocessing/test_file/"
    file_name = 'FlowDefinition_only_H.mff'
    
    # import mesh data from mff.file
    node_df, ele_df = MffProcessing.path_input(path + file_name)

    # aperture distribution
    trace_aperture = np.array([
        [[-5.70E-07, 2.25E-05, 9.78E-06], [-4.46E-07, 1.40E-05, 7.20E-07]]
    ])

    # 內寬分布
    aperture_dis = MeshApertureDistribution(node_df, ele_df, trace_aperture)
    # 已經內寬分布的element data frame
    cal_ele_df = aperture_dis.main_distribution()
    
    # head boundary 給node df 邊界條件
    # 生成DFN的 [X, Y, Z, SX, SY, SZ]
    region = np.array([10, 10, 10, 20, 20, 20])
    
    '''
    #   x
    # ax + by + cz + h (head)
    constant_10 = [0, 0, 0, 10]
    constant_0 = [0, 0, 0, 0]
    var_x_1 = [1, 0, 0, 0]
    boundary = np.array([
        constant_10,
        constant_0,
        var_x_1
    ])
    boundary_type = 'x'

    boundary_x_assign = HeadBoundaryAssign(node_df, boundary, boundary_type, region)
    # x 的水頭邊界
    cal_node_df_x = boundary_x_assign.main_assign(err=1e-3)

    #   y
    constant_10 = [0, 0, 0, 10]
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
    constant_10 = [0, 0, 0, 10]
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

    # 將以內寬分布df加水偷邊界node df 輸出
    # output mff.file of three boundaries
    # element
    MffProcessing.ele_output_mff(path + file_name, 'ele.mff', cal_ele_df)
    
    
    # node 
    MffProcessing.node_output_mff('ele.mff', 'liner_x_200.mff', cal_node_df_x)
    MffProcessing.node_output_mff('ele.mff', 'output_y_8.mff', cal_node_df_y)
    MffProcessing.node_output_mff('ele.mff', 'output_z_8.mff', cal_node_df_z)
    '''