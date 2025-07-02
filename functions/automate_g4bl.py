import subprocess
import os
import sys


# Set working directory to parent:
def setDirParent():
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
    sys.path.append(parent_dir)


# Create directory for next sim:
def createSimDir(dir):

    ''' Inputs:
    dir = name of new directory '''

    command = f"mkdir {dir}"
    process = subprocess.run(command, shell=True, capture_output=False)


# Copy nominal input file to sim directory:
def cpInput(file_path, dir):

    ''' Inputs:
    path = absolute path to template .in file
    dir = sim directory, where input file will be copied to '''

    command = f"cp {file_path} {dir}/g4bl_input_script.in"
    process = subprocess.run(command, shell=True, capture_output=False)


# Source bash profile:
def sourceBash(bash_path):

    ''' Inputs:
    bash_path = absolute path to .bash_profile '''

    command = f"source {bash_path}"
    process = subprocess.run(command, shell=True, capture_output=True, text=True)


# Run G4beamline:
def runG4bl(file_path):

    ''' Inputs:
    file_path = absolute path to g4bl input file
    
    Note that g4bl will write output to the directory in which the input file is located '''

    command = f"g4bl {file_path}"
    process = subprocess.run(command, shell=True, capture_output=True, text=True)
    return process.stdout, process.stderr


# Remove unnecessary output files:
def rmOutput(file_paths):

    ''' Inputs:
    file_paths = list of absolute paths to files to be removed '''

    for i in range(len(file_paths)):
        command = f"rm {file_paths[i]}"
        process = subprocess.run(command, shell=True, capture_output=False)
        del command, process