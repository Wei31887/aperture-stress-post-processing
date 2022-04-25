import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import time

from mesh_distribution import AperturePostProcessing, HeadBoundaryAssign, MeshApertureDistribution    
from mff_processing import MffProcessing
from stress_module import StressToAperture

if __name__ == '__main__':
    # import data file 
    # path = '/Users/weitim/Desktop/研究/研究主軸/FracMan_module/stress_on_aperture_distriution/'
    # file_name = 'test_x.mff'
    path = 'C:/Users/user/fracman_aperture_postprocessing/'
    file_name = 'test_x.mff'
    
    # import mesh data from mff.file
    node_df, ele_df = MffProcessing.path_input(path + file_name)

    # aperture distribution
    trace_aperture = np.array([
        [[-4 * (10**-7), 2 * (10**-5), 2 * (10**-5)], [-4 * (10**-7), 1 * (10**-5), 7 * (10**-7)]],
        [[-8 * (10**-7), 3 * (10**-5), 5 * (10**-7)], [-2 * (10**-7), 7 * (10**-6), 7 * (10**-7)]],
        [[-8 * (10**-7), 3 * (10**-5), 5 * (10**-7)], [-7 * (10**-7), 2 * (10**-5), 1 * (10**-5)]]
    ])
    
    aperture_dis = MeshApertureDistribution(node_df, ele_df, trace_aperture)
    cal_ele_df = aperture_dis.main_distribution()

    # head boundary
    region = np.array([5, 5, 5, 10, 10, 10])
    #   x
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
    cal_node_df_x = boundary_x_assign.main_assign()

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
    cal_node_df_y = boundary_y_assign.main_assign()

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
    cal_node_df_z = boundary_z_assign.main_assign()

    # output mff.file of three boundaries
    # element
    MffProcessing.ele_output_mff(path + file_name, 'ele.mff', cal_ele_df)
    
    # node 
    MffProcessing.node_output_mff(path + 'ele.mff', 'output_x.mff', cal_node_df_x)
    MffProcessing.node_output_mff(path + 'ele.mff', 'output_y.mff', cal_node_df_y)
    MffProcessing.node_output_mff(path + 'ele.mff', 'output_z.mff', cal_node_df_z)