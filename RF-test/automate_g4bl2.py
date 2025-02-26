import subprocess
import numpy as np

##### INPUTS #####

# Define directory where g4bl is installed:
g4bl_dir = '/Users/criggall/Documents/muon-cooling/' # <--- Location of .bash_profile

# Define working directory:
dir = '/Users/criggall/Documents/muon-cooling/RF-test/' # <-- All input files must be in this same directory

# Define file locations:
file = dir+'ref_finder_v1.in'
file_for_g4bl = '"'+file+'"'
det_file = dir+'detectors.txt'
beam_file = dir+'initial.dat'

# Specify whether to simulate a reference particle:
ref_particle = True

# Specify whether to simulate beam as well:
beam = True

# Specify beam type:
beam_type = 'gaussian' # <-- Options are 'gaussian' and 'file'

# Set constant parameters:
# ref_p = 225 # reference particle momentum (MeV/c)

# Define range for parameters to scan over:
# toffset = np.arange(0,3.01,0.01)
ref_p = np.arange(225,230.1,0.1)

# Set number of loops based on parameter scan space:
iterations = len(ref_p) # <-- Adding a second parameter to scan over will require a second loop below

##### FUNCTION DEFINITIONS #####

# Define function to write out G4bl output:
def write_out(out_dir):
    mkdir_command = f"mkdir {out_dir}"
    mkdir_process = subprocess.run(mkdir_command, shell=True, capture_output=False)

# Define function for modifying G4beamline input card:
def modify_g4bl_input(dir, file, beam_file, parameters, out_dir):

    with open(file, 'r') as f:
        lines = f.readlines()
        count = 0
        for i, line in enumerate(lines):
            
            # Adjust input file paths:
            if 'beam ascii' in line or 'beam gaussian' in line:
                if beam == True:
                    if beam_type == 'gaussian':
                        lines[i] = "beam gaussian particle=mu+ nEvents=10 beamZ=-700.0 sigmaX=100.0 sigmaY=100.0 sigmaXp=0.00 sigmaYp=0.00 meanMomentum=250.0 sigmaP=0.0 meanT=0.0 sigmaT=0\n"
                    if beam_type == 'file':
                        lines[i] = f"beam ascii file={beam_file} beamZ=$beamstart\n"
                # Replace with Gaussian beam to reduce sim time if False:
                else:
                    lines[i] = "beam gaussian particle=mu+ nEvents=100 beamZ=-700.0 sigmaX=10.0 sigmaY=10.0 sigmaXp=0.100 sigmaYp=0.100 meanMomentum=200.0 sigmaP=4.0 meanT=0.0 sigmaT=0.0\n"
            elif 'include' in line and 'detectors.txt' in line:
                lines[i] = f"include {dir}detectors.txt\n"
            elif 'include' in line and 'solangles1.dat' in line:
                lines[i] = f"include {dir}solangles1.dat\n"
            elif 'include' in line and 'solangles2.dat' in line:
                lines[i] = f"include {dir}solangles2.dat\n"
            elif 'include' in line and 'abs_place7_31.txt' in line:
                lines[i] = f"include {dir}abs_place7_31.txt\n"
            elif 'include' in line and 'RFplace7_31.txt' in line:
                lines[i] = f"include {dir}RFplace7_31.txt\n"
            elif 'include' in line and 'sol_place7_31.txt' in line:
                lines[i] = f"include {dir}sol_place7_31.txt\n"
            elif 'virtualdetector' in line:
                lines[i] = f"virtualdetector DetLast file={out_dir}outlast format=ascii radius=300 color=0,1,0 length=0.001 material=Vacuum\n"

            # Check if reference particle is active and remove if False:
            if 'trace' in line:
                count += 1
                if ref_particle == False:
                    lines[i] = "\n"
            if 'reference' in line:
                if ref_particle == False:
                    lines[i] = "\n"
            
        # Add reference particle if True:
        if ref_particle == True and count == 0:
            lines.append(f'reference referenceMomentum={ref_p} particle=mu+ beamZ=0.0\n')
            lines.append(f'trace nTrace=1 format=ascii file="{g4bl_dir}TraceParticle.txt"\n')
        
    with open(file, 'w') as f:
        f.writelines(lines)

# Define function to modify detector file:
def modify_detector_file(out_dir, det_file):
        
        with open(det_file, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                
                # Adjust output file paths:
                if 'file=' in line:
                    lines[i] = f"virtualdetector Det{i+1} file={out_dir}out{i+1} format=ascii radius=360 color=0,1,0 length=0.001 material=Vacuum\n"
        
        with open(det_file, 'w') as f:
            f.writelines(lines)

# Define function to run G4beamline:
def run_g4bl_sim(file):
    g4bl_command = f"source /Users/criggall/Documents/muon-cooling/.bash_profile && g4bl {file_for_g4bl}"
    g4bl_process = subprocess.run(g4bl_command, shell=True, capture_output=True, text=True)
    return g4bl_process.stdout, g4bl_process.stderr

# Define function to copy .in card to output directory:
def copy_in(file, out_dir):
    copy_in_command = f"cp {file} {out_dir}/G4blInputCard.in"
    copy_in_process = subprocess.run(copy_in_command, shell=True, capture_output=False)

# Define function to move reference particle output and delete additional trace files:
def mv_ref_file(g4bl_dir, out_dir):
    mv_ref_file_command = f"mv {g4bl_dir}ReferenceParticle.txt {out_dir} && rm {g4bl_dir}TuneParticle.txt && rm {g4bl_dir}TraceParticle.txt"
    mv_ref_file_process = subprocess.run(mv_ref_file_command, shell=True, capture_output=False)

# Define function to remove detector output files for reference-particle-only sim:
def rm_det_out(out_dir):
    rm_det_out_command = f"rm {out_dir}out*.txt"
    rm_det_out_process = subprocess.run(rm_det_out_command, shell=True, capture_output=False)

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
    
    # Copy .in card to output directory:
    copy_in(file, out_dir)

    # Set configurable parameters:
    parameters = {
        'ref_p' : ref_p[j]
    }

    # Execute functions:
    modify_g4bl_input(dir, file, beam_file, parameters, out_dir)
    if ref_particle == False:
        modify_detector_file(out_dir, det_file)
    stdout, stderr = run_g4bl_sim(file)
    if ref_particle == True:
        mv_ref_file(g4bl_dir, out_dir)
    if beam == False:
        rm_det_out(out_dir)
        print(f'Removing detector files from directory {out_dir}')
    
    if stderr:
            print("Error running simulation:", stderr)
    else:
            print("Simulation completed successfully")
            print("Output:", stdout)