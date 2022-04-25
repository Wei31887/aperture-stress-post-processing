import numpy as numpy
import pandas as pd

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
            c_idx = 0 
            for row in range(len(lines)):
                if row in range(start_row_idx+1, end_row_idx):
                    f_out.write(
                        f'{change_lines[c_idx][0]:>8s}'
                        f'{change_lines[c_idx][1]:>8s}'
                        f'{change_lines[c_idx][2]:>8s}'
                        f'{change_lines[c_idx][3]:>8s}'
                        f'{change_lines[c_idx][4]:>8s}'
                        f'{change_lines[c_idx][5]:>6s}'
                        f'{change_lines[c_idx][6]:>13s}'
                        f'{change_lines[c_idx][7]:>13s}'
                        f'{change_lines[c_idx][8]:>13s}\n'
                    )
                    c_idx += 1 
                else:  
                    f_out.write(orin_lines[row])    
               
    @staticmethod
    def node_output_mff(data_path, file_name, cal_node_df):
        with open(file_name, mode='w') as f_out, open(data_path, mode='r') as f_copy:
            copy_content = f_copy.readlines()

            # get change part: node
            for idx, line in enumerate(copy_content):
                if 'NODE' in line:
                    start_row_idx = idx+1
                if 'FracElem' in line:
                    end_row_idx = idx-1

            # change lines with rewrite node data
            change_lines = []
            for df_idx, row_idx in enumerate(range(start_row_idx, end_row_idx)):
                change_line = copy_content[row_idx]
                change_line = change_line.strip()
                change_line = change_line.split()
                # change part: type, head, grp
                #   type
                change_line[4] = f"{int(cal_node_df['Type'].iloc[df_idx])}"
                #   head
                change_line[5] = f"{cal_node_df['H'].iloc[df_idx]:.4E}"
                #   Grp
                change_line[7] = f"{int(cal_node_df['Grp'].iloc[df_idx])}"
                change_lines.append(change_line)

            # write the change part
            change_row_idx = 0
            for row_idx in range(len(copy_content)):
                if row_idx in range(start_row_idx, end_row_idx):
                    f_out.write(
                        f'{change_lines[change_row_idx][0]:>8s}'
                        f'{change_lines[change_row_idx][1]:>19s}'
                        f'{change_lines[change_row_idx][2]:>19s}'
                        f'{change_lines[change_row_idx][3]:>19s}'
                        f'{change_lines[change_row_idx][4]:>3s}'
                        f'{change_lines[change_row_idx][5]:>12s}'
                        f'{change_lines[change_row_idx][6]:>12s}'
                        f'{change_lines[change_row_idx][7]:>9s}\n'
                    )
                    change_row_idx += 1
                else:
                    f_out.write(copy_content[row_idx])

    @staticmethod
    def LR_CRLF_changer(file_name):
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
                    

if __name__ == '__main__':
    path = 'C:/Users/user/fracman_aperture_postprocessing/'
    file_name = 'test_y.mff'
    
    # import mesh data from mff.file
    node_df, ele_df = MffProcessing.path_input(path + file_name)
    
    MffProcessing.node_output_mff(path + file_name, 'x.mff', node_df)