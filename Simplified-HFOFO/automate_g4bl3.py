import subprocess
import numpy as np

##### INPUTS #####

# Define directory where g4bl is installed:
g4bl_dir = '/Users/criggall/Documents/muon-cooling/' # <--- Location of .bash_profile

# Define working directory:
dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/' # <-- All input files must be in this same directory

# Define file locations:
file = dir+'g4bl_input_card.in'
file_for_g4bl = '"'+file+'"'

# Define range of solenoid current scan:
sol_curr = np.arange(90,100,0.1)

# Set number of loops based on scan space:
iterations = len(sol_curr)

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
            
            # Adjust solenoid constant current:
            if 'param sol_current' in line:
                lines[i] = f"param sol_current={parameters['sol_curr']} #amps/mm^2\n"
        
    with open(file, 'w') as f:
        f.writelines(lines)

# Define function to run G4beamline:
def run_g4bl_sim(file):
    g4bl_command = f"source /Users/criggall/Documents/muon-cooling/.bash_profile && g4bl {file_for_g4bl}"
    g4bl_process = subprocess.run(g4bl_command, shell=True, capture_output=True, text=True)
    return g4bl_process.stdout, g4bl_process.stderr

# Define function to copy .in card to output directory:
def copy_in(file, out_dir):
    copy_in_command = f"cp {file} {out_dir}/g4bl_input_card.in"
    copy_in_process = subprocess.run(copy_in_command, shell=True, capture_output=False)

# Define function to move reference particle output and delete additional trace files:
def mv_ref_file(g4bl_dir, out_dir):
    mv_ref_file_command = f"mv ReferenceParticle.txt {out_dir} && rm {g4bl_dir}TuneParticle.txt && rm {g4bl_dir}coil.dat"
    mv_ref_file_process = subprocess.run(mv_ref_file_command, shell=True, capture_output=False)

##### MAIN LOOP #####

for j in range(iterations):

    print(f"Now running simulation {j+1}...")

    # Define output directory:
    out_dir = dir+f'g4bl-output-sim{j+1}/'

    # Make output directory if it does not already exist:
    dir_exists_command = f'if test -d {out_dir}; then echo 1; fi'
    dir_exists_process = subprocess.run(dir_exists_command, shell=True, capture_output=True)
    dir_exists_output = str(dir_exists_process.stdout)[2]
    if dir_exists_output != '1':
        write_out(out_dir)

    # Set configurable parameters:
    parameters = {
        'sol_curr' : sol_curr[j]
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