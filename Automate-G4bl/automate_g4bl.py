import subprocess

# .in file location:
dir = '/Users/criggall/Documents/muon-cooling/Automate-G4bl/' # <-- All input files must be in this same directory!
file = dir+'track_v7.in'
file_for_g4bl = '"'+file+'"'
det_file = dir+'detectors.txt'
beam_file = dir+'initial.dat'

# Set configurable parameters:
parameters = {
    'beamstart':-700, # adjusts initial z position offset
    'beamtime':-0.671, # adjusts initial time offset
    'density':0.014 # GH2 density
}

# Define function for modifying G4beamline input card:
def modify_g4bl_input(dir, file, beam_file, parameters):

    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            
            # Adjust input file paths:
            if 'beam ascii file=' in line:
                 lines[i] = f"beam ascii file={beam_file} beamZ=$beamstart\n"
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

            # Update configurable paramters:
            if 'param beamstart' in line:
                lines[i] = f"param beamstart={parameters['beamstart']}\n"
            elif 'param beamtime' in line:
                lines[i] = f"param beamtime={parameters['beamtime']}\n"
            elif 'material GH2' in line:
                lines[i] = f"material GH2 Z=1 A=1.01 density={parameters['density']}\n"
        
    with open(file, 'w') as f:
        f.writelines(lines)

def modify_detector_file(dir, det_file):
    
    with open(det_file, 'r') as f:
         lines = f.readlines()
         for i, line in enumerate(lines):
              
              # Adjust output file paths:
              if 'file=' in line:
                   lines[i] = f"virtualdetector Det{i+1} file={dir}out{i+1} format=ascii radius=360 color=0,1,0 length=0.001 material=Vacuum\n"
    
    with open(det_file, 'w') as f:
        f.writelines(lines)

# Define function to run G4beamline:
def run_g4bl_sim(file):
    g4bl_command = f"source .bash_profile && g4bl {file_for_g4bl}"
    g4bl_process = subprocess.run(g4bl_command, shell=True, capture_output=True, text=True)
    return g4bl_process.stdout, g4bl_process.stderr

# Execute functions:
modify_g4bl_input(dir, file, beam_file, parameters)
modify_detector_file(dir, det_file)
stdout, stderr = run_g4bl_sim(file)

if stderr:
        print("Error running simulation:", stderr)
else:
        print("Simulation completed successfully")
        print("Output:", stdout)