import math

# Header file for block processing, includes parameters and paths
isCluster = True

folder_path = data_path + sample_name + "/" + outp_ID + "/"

# folder to which the current code is copied and where it is run from
code_run_path = folder_path + "code/"
slurm_path = folder_path + "slurm_files/"
# these are created by the preparations folder
error_path = folder_path+"error_files/"
output_path = folder_path+"output_files/"

# create file no save n_comp (number of cc3d components) for each procesed block
n_comp_filepath         = folder_path+"n_comp.txt"
step01_info_filepath  = folder_path+"step01_info.txt"
step02A_info_filepath = folder_path+"step02A_info.txt"
step02B_info_filepath = folder_path+"step02B_info.txt"
step03_info_filepath  = folder_path+"step03_info.txt"

# memory need per block (in MB)
memory_step01_number = int(1.1*max_bs_z*max_bs_y*max_bs_x*(8+8+8)/1000/1000)
memory_step01 = str(memory_step01_number)
memory_step02A = str(int(memory_step01_number*0.05))
memory_step02B = str(int(memory_step01_number*0.5))
memory_step03 = str(int(memory_step01_number*1))

run_hours = str(2)

max_labels_block = max_bs_z*max_bs_y*max_bs_x
