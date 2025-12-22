import subprocess
import os
import sys
import numpy as np

########## FUNCTIONS ##########

# Define function to write out G4bl output:
def write_out(out_dir):
    mkdir_command = f"mkdir {out_dir}"
    mkdir_process = subprocess.run(mkdir_command, shell=True, capture_output=False)

# Define function to copy input card to output directory:
def copy_in(file, out_dir):
    copy_in_command = f"cp {file} {out_dir}/g4bl_input_script.in"
    copy_in_process = subprocess.run(copy_in_command, shell=True, capture_output=False)

# Define function to run G4beamline:
def run_g4bl_sim(file):
    g4bl_command = f"g4bl {file_for_g4bl}"
    g4bl_process = subprocess.run(g4bl_command, shell=True, capture_output=True, text=True)
    return g4bl_process.stdout, g4bl_process.stderr

# Define function to move output files to simulation directory:
def mv_ref_file(dir, out_dir):
    #mv_ref_file_command = f"mv {dir}*.txt {out_dir} & mv g4bl.out {out_dir} & mv {dir}*.dat"
    # Update this command to use the specific names of your output files to avoid accidentally moving other files!
    mv_ref_file_command = f"mv AllTracks.txt {out_dir}"
    mv_ref_file_process = subprocess.run(mv_ref_file_command, shell=True, capture_output=False)

########## DEFINE SCAN PARAMETERS ##########

# Define file locations:
file = 'track_v7.in'

file_for_g4bl = '"'+file+'"'

# Define range for parameters to scan over:
# ref_p = np.arange(230,250,1)
# ref_p = np.arange(237,246,0.1)
# ref_p = np.arange(245,246,0.01)
# ref_p = np.arange(245.50,247.50,0.01)
ref_p = np.arange(246.38,246.40,0.001)

# Set number of loops based on parameter scan space:
iterations = len(ref_p)

# Define function for modifying G4beamline input card:
def modify_g4bl_input(dir, file, parameters, out_dir):

    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # Update configurable parameters:
            if 'param ref_p=' in line:
               lines[i] = f"param ref_p={parameters['reference momentum']}\n"
        
    with open(file, 'w') as f:
        f.writelines(lines)

##### MAIN LOOP #####

for j in range(iterations):

    print(f"Now running simulation {j+1}...")

    # Define output directory:
    out_dir = f'sim{j+1}/'

    # Make output directory if it does not already exist:
    dir_exists_command = f'if test -d {out_dir}; then echo 1; fi'
    dir_exists_process = subprocess.run(dir_exists_command, shell=True, capture_output=True)
    dir_exists_output = str(dir_exists_process.stdout)[2]
    if dir_exists_output != '1':
        write_out(out_dir)

    # Set configurable parameters:
    parameters = {
        'reference momentum' : ref_p[j]
    }

    # Modify input files:
    modify_g4bl_input(dir, file, parameters, out_dir)

    # Copy input card to simulation directory:
    copy_in(file, out_dir)

    # Execute g4bl:
    stdout, stderr = run_g4bl_sim(file)
    
    # Print output to terminal:
    if stderr:
            print("Error running simulation:", stderr)
    else:
            print("Simulation completed successfully")
            print("Output:", stdout)
    
    # Move output files to simulation directory:
    mv_ref_file(dir, out_dir)