import numpy as np
import pandas as pd
import math 
import shapely.geometry
import shapely.ops
import geopandas as gpd
import matplotlib.pyplot as plt
import random

# mff. inputer & outputer
class MffProcessing(object):  
    @staticmethod
    def path_input(data_path):
        """ Import .mff file from FracMan and get the dataframe of node & element 
        """
        with open(data_path, mode="r") as f:
            lines = f.readlines()
            lines = [l.strip() for l in lines]

            data_start_row = -1

            # Block1: get rid of HEADER
            for row_idx, line in enumerate(lines):
                if "***" in line:
                    data_start_row = row_idx + 2
                    break
            
            columns1 = lines[data_start_row-1].split(" ")
            columns1 = [c for c in columns1 if c]

            # Block2: get node data
            data_node = []
            node_start_row = data_start_row
            for row_idx in range(data_start_row, len(lines)):
                text_list = lines[row_idx].split(" ")
                try:
                    text_list = [float(t) for t in text_list if t]
                    data_node.append(text_list)
                except Exception as e:
                    ele_start_row  = row_idx
                    break

            # Block3: get element data
            columns2 = lines[ele_start_row].split(' ')
            columns2 = [c for c in columns2[:-1] if c]
            
            data_ele = []
            for row_idx in range(ele_start_row + 1, len(lines)):
                text_list = lines[row_idx].split(' ')
                try:
                    text_list = [float(t) for t in text_list if t]
                    data_ele.append(text_list)
                except Exception:
                    break

            data_node_df = pd.DataFrame(data_node, columns=columns1)
            data_node_df.set_index('NODE', inplace=True)
            data_ele_df = pd.DataFrame(data_ele, columns=columns2)
            data_ele_df.set_index('FracElem', inplace=True)
            
            return data_node_df, data_ele_df

    @staticmethod
    def ele_output_mff(data_path, file_name, cal_ele_df):
        with open(file_name, mode='w') as f_out, \
            open(data_path, mode='r') as f_copy:
            
            orin_lines = f_copy.readlines()
            lines = [l.strip() for l in orin_lines]

            # Get change part index: element part
            start_row_idx = 0
            for idx, line in enumerate(lines):
                if 'FracElem' in line:
                    start_row_idx = idx
                if 'MatrixElem' in line:
                    end_row_idx = idx

            # Change the data with calculated one
            change_lines = []
            data_idx = 0
            for change_row in range(start_row_idx + 1, end_row_idx):
                change_line = lines[change_row].strip()
                change_line = change_line.split(' ')
                change_line = [c for c in change_line if c]
                state = True
                for columns_idx, data in enumerate(change_line):
                    if columns_idx == 6 and state:
                        change_line[columns_idx] = f"{cal_ele_df['Trans'].iloc[data_idx]:.4E}"
                        state = False
                    if columns_idx == 8:
                        change_line[columns_idx] = f"{cal_ele_df['Apert'].iloc[data_idx]:.4E}"
                data_idx += 1    
                change_lines.append(change_line)

            # Write the output file
            out_lines = []
            c_idx = 0 
            for row in range(len(lines)):
                if row in range(start_row_idx+1, end_row_idx):
                    f_out.write(f'{change_lines[c_idx][0]:>8s}'
                                f'{change_lines[c_idx][1]:>8s}'
                                f'{change_lines[c_idx][2]:>8s}'
                                f'{change_lines[c_idx][3]:>8s}'
                                f'{change_lines[c_idx][4]:>8s}'
                                f'{change_lines[c_idx][5]:>6s}'
                                f'{change_lines[c_idx][6]:>13s}'
                                f'{change_lines[c_idx][7]:>13s}'
                                f'{change_lines[c_idx][8]:>13s}\n')
                    c_idx += 1 
                else:  
                    f_out.write(orin_lines[row])  
        
        # Turn into CRLF formate
        # replacement strings
        WINDOWS_LINE_ENDING = b'\r\n'
        UNIX_LINE_ENDING = b'\n'

        with open(file_name, 'rb') as open_file:
            content = open_file.read()

        # Unix ➡ Windows
        content = content.replace(UNIX_LINE_ENDING, WINDOWS_LINE_ENDING)

        with open(file_name, 'wb') as open_file:
            open_file.write(content)     
               
    @staticmethod
    def node_output_mff(data_path, file_name, cal_node_df):
        with open(file_name, mode='w') as f_out, open(data_path, mode='r') as f_copy:
            copy_content = f_copy.readlines()
                            
class HeadBoundaryAssign(object):
    def __init__(self, node_df, boundary, boundary_type, region):
        self.node_df = node_df
        self.boundary = boundary
        self.boundary_type = boundary_type
        self.region = region
        self.set_boundary()
        
    def set_boundary(self):
        if self.boundary_type == 'x':
            self.boundary_for_set = [
                ([self.region[0]+self.region[3]/2, 0, 0], self.boundary[0]),
                ([self.region[0]-self.region[3]/2, 0, 0], self.boundary[1]),
                ([0, self.region[1]+self.region[4]/2, 0], self.boundary[2]),
                ([0, self.region[1]-self.region[4]/2, 0], self.boundary[2]),
                ([0, 0, self.region[2]+self.region[5]/2], self.boundary[2]),
                ([0, 0, self.region[2]-self.region[5]/2], self.boundary[2]),
            ]
        elif self.boundary_type == 'y':
            self.boundary_for_set = [
                ([self.region[0]+self.region[3]/2, 0, 0], self.boundary[2]),
                ([self.region[0]-self.region[3]/2, 0, 0], self.boundary[2]),
                ([0, self.region[1]+self.region[4]/2, 0], self.boundary[0]),
                ([0, self.region[1]-self.region[4]/2, 0], self.boundary[1]),
                ([0, 0, self.region[2]+self.region[5]/2], self.boundary[2]),
                ([0, 0, self.region[2]-self.region[5]/2], self.boundary[2]),
            ]
        elif self.boundary_type == 'z':
            self.boundary_for_set = [
                ([self.region[0]+self.region[3]/2, 0, 0], self.boundary[2]),
                ([self.region[0]-self.region[3]/2, 0, 0], self.boundary[2]),
                ([0, self.region[1]+self.region[4]/2, 0], self.boundary[2]),
                ([0, self.region[1]-self.region[4]/2, 0], self.boundary[2]),
                ([0, 0, self.region[2]+self.region[5]/2], self.boundary[0]),
                ([0, 0, self.region[2]-self.region[5]/2], self.boundary[1]),
            ]
            
    def function_give_head(self, x, y, z, boundary_type):
        ''' Head = hx*x + hy*y + hz*z + h0'''
        hx = boundary_type[1][0]
        hy = boundary_type[1][1]
        hz = boundary_type[1][2]
        h0 = boundary_type[1][3]
        head = x*hx + y*hy + z*hz + h0

        return head
    
    def main_assign(self, err=1e-8):
        '''
        Set the node data the given head boundary.
        read through the node data and determine the node is on the boundary or not
        '''    
        output_node_df = self.node_df.copy() 
        # read through the node df and give the corespond head
        for node in output_node_df.index:
            x = output_node_df['X'].loc[node]
            y = output_node_df['Y'].loc[node]
            z = output_node_df['Z'].loc[node]
            # determine the node is on boundary or not
            # boundary1
            if abs(x - self.boundary_for_set[0][0][0]) <= err:
                tem_head = self.function_give_head(x, y, z, self.boundary_for_set[0])
                output_node_df['H'].loc[node] = tem_head
                output_node_df['Type'].loc[node] = 1
                if self.boundary_type == 'x':
                    output_node_df['Grp'].loc[node] = 1
                elif self.boundary_type == 'y':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'z':
                    output_node_df['Grp'].loc[node] = 3
            # boundary2
            elif abs(x - self.boundary_for_set[1][0][0]) <= err:
                tem_head = self.function_give_head(x, y, z, self.boundary_for_set[1])
                output_node_df['H'].loc[node] = tem_head
                output_node_df['Type'].loc[node] = 1
                if self.boundary_type == 'x':
                    output_node_df['Grp'].loc[node] = 2
                elif self.boundary_type == 'y':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'z':
                    output_node_df['Grp'].loc[node] = 3
            # boundary3
            elif abs(y - self.boundary_for_set[2][0][1]) <= err:
                tem_head = self.function_give_head(x, y, z, self.boundary_for_set[2])
                output_node_df['H'].loc[node] = tem_head
                output_node_df['Type'].loc[node] = 1
                if self.boundary_type == 'x':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'y':
                    output_node_df['Grp'].loc[node] = 1
                elif self.boundary_type == 'z':
                    output_node_df['Grp'].loc[node] = 3
            # boundary4
            elif abs(y - self.boundary_for_set[3][0][1]) <= err:
                tem_head = self.function_give_head(x, y, z, self.boundary_for_set[3])
                output_node_df['H'].loc[node] = tem_head
                output_node_df['Type'].loc[node] = 1
                if self.boundary_type == 'x':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'y':
                    output_node_df['Grp'].loc[node] = 2
                elif self.boundary_type == 'z':
                    output_node_df['Grp'].loc[node] = 3
            # boundary5
            elif abs(z - self.boundary_for_set[4][0][2]) <= err:
                tem_head = self.function_give_head(x, y, z, self.boundary_for_set[4])
                output_node_df['H'].loc[node] = tem_head
                output_node_df['Type'].loc[node] = 1
                if self.boundary_type == 'x':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'y':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'z':
                    output_node_df['Grp'].loc[node] = 1
            # boundary6
            elif abs(z - self.boundary_for_set[5][0][2]) <= err:
                tem_head = self.function_give_head(x, y, z, self.boundary_for_set[5])
                output_node_df['H'].loc[node] = tem_head
                output_node_df['Type'].loc[node] = 1
                if self.boundary_type == 'x':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'y':
                    output_node_df['Grp'].loc[node] = 3
                elif self.boundary_type == 'z':
                    output_node_df['Grp'].loc[node] = 2
            else:
                output_node_df['Type'].loc[node] = 0
                output_node_df['Grp'].loc[node] = 0
                pass     
        
        return output_node_df

class MeshApertureDistribution(object):
    def __init__(self, node_df, ele_df, aperture_func):
        self.node_df = node_df
        self.ele_df = ele_df 
        self.aperture_func = aperture_func
    
    def merge_meshes(self, node_list):
        """ Merge the meshes of one single plane into one polygon 
        """
        fracture_polygon = []
        for nodes in node_list:
            tem_polygon = shapely.geometry.Polygon([
                    (self.node_df['X'][nodes[0]], self.node_df['Y'][nodes[0]], self.node_df['Z'][nodes[0]]),
                    (self.node_df['X'][nodes[1]], self.node_df['Y'][nodes[1]], self.node_df['Z'][nodes[1]]),
                    (self.node_df['X'][nodes[2]], self.node_df['Y'][nodes[2]], self.node_df['Z'][nodes[2]])
                    ])
            # fracture_polygon.append(tem_polygon)
            merged_polygon = shapely.ops.unary_union([tem_polygon])

        return merged_polygon
    
    def collect_node_data(self, fracture_index):    
        # Collect mesh data from one fracture
        element_list = [i for i in self.ele_df['Set#'].index \
                            if self.ele_df['Set#'][i] == fracture_index]
        node_list = []
        for element in element_list:
            node_list.append((
                self.ele_df['Node1'][element],
                self.ele_df['Node2'][element],
                self.ele_df['Node3'][element]
                             ))
        return node_list
    
    def polygon_tracelength(self, polygon):
        """ Determine the trace length of one fracture polygon
        """
        polygon_exterior = [exterior for exterior in polygon.exterior.coords]
        cal_exterior = polygon_exterior[1:]
        length_list = []
        for coords in polygon_exterior:
            for cal_coords in cal_exterior:
                length = ((coords[0] - cal_coords[0])**2 \
                            + (coords[1] - cal_coords[1])**2 \
                            + (coords[2] - cal_coords[2])**2) ** 0.5
                length_list.append(length)
            if len(cal_exterior) == 1:
                break
            cal_exterior = cal_exterior[1:]
        two_times_radius = max(length_list)
        # print(len(length_list))
        trace_length = (math.pi**0.5) * (two_times_radius / 2)
        return trace_length
    
    def aperture_distribution(self, fracture_index, trace_length):
        mean_func = self.aperture_func[fracture_index-1][0]
        std_func = self.aperture_func[fracture_index-1][1]
        aperture_mean = (mean_func[0]**2) * trace_length + mean_func[1]*trace_length + mean_func[2]
        aperture_std = (std_func[0]**2) * trace_length + std_func[1]*trace_length + std_func[2]
        
        return aperture_mean, aperture_std
    
    def main_distribution(self):
        """ 
        1. Read though the df of element
            1.1. Find the mesh belong to one single plane
            1.2. Merge into one polygon
            1.3. Determine the trace length of plane
            1.4. Apply the corresponding aperture
        """
        # Read though the df of element
        fracture_index = 1 
        fracture_set_max = max(self.ele_df['Set#'])
        output_ele_df = self.ele_df
        while True:
            # 1.1. Find the mesh belong to one single plane
            tem_node_list = self.collect_node_data(fracture_index)
            
            # 1.2. Merge into one polygon
            merge_polygon = self.merge_meshes(tem_node_list)
            # print(merge_polygon)
        
            # 1.3. Determine the trace length of plane
            trace_length = self.polygon_tracelength(merge_polygon)
            
            # 1.4. Apply the corresponding aperture
            # And the corresponding trans, stor
            aperture_mean, aperture_std = self.aperture_distribution(fracture_index, trace_length)
            for ele in self.ele_df[self.ele_df['Set#'] == fracture_index].index:
                tem_aperture = random.gauss(aperture_mean, aperture_std)
                output_ele_df['Apert'][ele] = tem_aperture
                output_ele_df['Trans'][ele] = 817500 * tem_aperture**3
            fracture_index += 1
            if fracture_index > fracture_set_max:
                return output_ele_df

class AperturePostProcessing(HeadBoundaryAssign, MeshApertureDistribution):
    def __init__(self, node_df, ele_df, aperture_func, boundary, boundary_type, region):
        HeadBoundaryAssign.__init__(self, node_df, boundary, boundary_type, region)
        MeshApertureDistribution.__init__(self, node_df, ele_df, aperture_func)

if __name__ == '__main__':
    pass
