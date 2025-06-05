import subprocess
import numpy as np

##### INPUTS #####

# Define directory where g4bl is installed:
g4bl_dir = '/Users/criggall/Documents/muon-cooling/Solenoid-Study/' # Location of .bash_profile

# Define working directory:
dir = g4bl_dir+'single-coil/'

# Define file location:
# file = dir+'singlecoil.in'
file = dir+'secondcoil.in'
file_for_g4bl = '"'+file+'"'

# # Define range of coil length scan:
# coil_lengths = np.arange(5,400,5)

# Define range of coil spacing scan:
periods = np.arange(100,1500,100)

# Set number of loops based on scan space:
# iterations = len(coil_lengths)
iterations = len(periods)

##### FUNCTION DEFINITIONS #####

# Define function to write out G4bl output:
def write_out(out_dir):
    mkdir_command = f"mkdir {out_dir}"
    mkdir_process = subprocess.run(mkdir_command, shell=True, capture_output=False)

# Define function for modifying G4beamline input card:
def modify_g4bl_input(dir, file, parameters, out_dir):

    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):

            # if 'param sollength=' in line:
            #     lines[i] = f"param sollength={parameters['coil length']}\n"
            if 'param period=' in line:
                lines[i] = f"param period={parameters['period']}\n"
        
    with open(file, 'w') as f:
        f.writelines(lines)

# Define function to run G4beamline:
def run_g4bl_sim(file):
    g4bl_command = f"source {g4bl_dir}.bash_profile && g4bl {file_for_g4bl}"
    g4bl_process = subprocess.run(g4bl_command, shell=True, capture_output=True, text=True)
    return g4bl_process.stdout, g4bl_process.stderr

# Define function to copy .in card to output directory:
def copy_in(file, out_dir):
    copy_in_command = f"cp {file} {out_dir}/g4bl_input_card.in"
    copy_in_process = subprocess.run(copy_in_command, shell=True, capture_output=False)

# Define function to move reference particle output and delete additional trace files:
def mv_ref_file(g4bl_dir, out_dir):
    mv_ref_file_command = f"mv AllTracks.txt {out_dir} && rm det.txt && rm Coil1.dat"
    mv_ref_file_process = subprocess.run(mv_ref_file_command, shell=True, capture_output=False)

##### MAIN LOOP #####

for j in range(iterations):

    print(f"Now running simulation {j+1}...")

    # Define output directory:
    # out_dir = dir+f'coil_length_scan/g4bl-output-sim{j+1}/'
    out_dir = dir+f'coil_spacing_coarse_scan/g4bl-output-sim{j+1}/'
    # out_dir = dir+f'coil_spacing_fine_scan/g4bl-output-sim{j+1}/'

    # Make output directory if it does not already exist:
    dir_exists_command = f'if test -d {out_dir}; then echo 1; fi'
    dir_exists_process = subprocess.run(dir_exists_command, shell=True, capture_output=True)
    dir_exists_output = str(dir_exists_process.stdout)[2]
    if dir_exists_output != '1':
        write_out(out_dir)

    # Set configurable parameters:
    parameters = {
        'period' : periods[j]
    }

    # Modify input files:
    modify_g4bl_input(dir, file, parameters, out_dir)

    # Copy .in card to output directory:
    copy_in(file, out_dir)

    # Execute g4bl:
    stdout, stderr = run_g4bl_sim(file)
    
    # Organize output files:
    mv_ref_file(g4bl_dir, out_dir)
    
    # Print output to terminal:
    if stderr:
            print("Error running simulation:", stderr)
    else:
            print("Simulation completed successfully")
            print("Output:", stdout)