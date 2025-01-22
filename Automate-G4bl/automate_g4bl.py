import subprocess

# .in file location:
file = '/Users/criggall/Documents/muon-cooling/Automate-G4bl/track_v7.in'
file_for_g4bl = '"'+file+'"'

# Set configurable parameters:
parameters = {
    'beamstart':-500, # adjusts initial z position offset
    'beamtime':-0.671, # adjusts initial time offset
    'density':0.014 # GH2 density
}

# Define function for modifying G4beamline input card:
def modify_g4bl_input(file, parameters):

    with open(file, 'r') as f:
        lines = f.readlines()

        for i, line in enumerate(lines):
            if 'param beamstart' in line:
                lines[i] = f"param beamstart={parameters['beamstart']}\n"
            elif 'param beamtime' in line:
                lines[i] = f"param beamtime={parameters['beamtime']}\n"
            elif 'gasdensity' in line:
                lines[i] = f"material GH2 Z=1 A=1.01 density={parameters['density']}\n"
        
    with open(file, 'w') as f:
        f.writelines(lines)

# Define function to run G4beamline:
def run_g4bl_sim(file):
    g4bl_command = f"source .bash_profile && g4bl {file_for_g4bl}"
    g4bl_process = subprocess.run(g4bl_command, shell=True, capture_output=True, text=True)
    return g4bl_process.stdout, g4bl_process.stderr

# Execute functions:
modify_g4bl_input(file, parameters)
stdout, stderr = run_g4bl_sim(file)

if stderr:
        print("Error running simulation:", stderr)
else:
        print("Simulation completed successfully")
        print("Output:", stdout)