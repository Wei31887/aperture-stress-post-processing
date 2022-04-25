import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

from mesh_distribution import MffProcessing, AperturePostProcessing    
from stress_module import StressToAperture

if __name__ == '__main__':
    path = '/Users/weitim/Desktop/研究/研究主軸/FracMan_module/stress_on_aperture_distriution/'
    file_name = 'test_x.mff'
    
    # import mesh data from mff.file
    node_df, ele_df = MffProcessing.path_input(path + file_name)

    # aperture distribution
    trace_aperture = np.array([
        [[-4 * (10**-7), 2 * (10**-5), 2 * (10**-5)], [-4 * (10**-7), 1 * (10**-5), 7 * (10**-7)]],
        [[-8 * (10**-7), 3 * (10**-5), 5 * (10**-7)], [-2 * (10**-7), 7 * (10**-6), 7 * (10**-7)]],
        [[-8 * (10**-7), 3 * (10**-5), 5 * (10**-7)], [-7 * (10**-7), 2 * (10**-5), 1 * (10**-5)]]
    ])
    
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
    
    mesh_dis = AperturePostProcessing(node_df, ele_df, trace_aperture, boundary, boundary_type, region)
    cal_ele_df = mesh_dis.main_distribution()
    cal_node_df = mesh_dis.main_assign()
    
    print(cal_ele_df.head(5))
    print(cal_node_df.head(5))