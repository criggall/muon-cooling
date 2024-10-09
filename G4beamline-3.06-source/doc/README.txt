                             G4beamline 3.06
                              by Tom Roberts
               Copyright (C) 2003-2018 by Tom Roberts, Muons, Inc.
                             All rights reserved.

                       http://g4beamline.muonsinc.com


LICENSE
-------

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

    http://www.gnu.org/copyleft/gpl.html

NOTE: This program uses several open-source libraries. Their license terms
can be found in the Acknowledgements section of G4beamlineUsersGuide.pdf 
(in the doc directory of the distribution).


GENERAL
-------

The general reference for G4beamline is the User's Guide, located in the
doc directory named G4beamlineUsersGuide.pdf.

It is also available on the web: http://g4beamline.muonsinc.com


INSTALLATION
------------

There is a "Quick Start" section in the User's Guide.


INSTALLING THE REQUIRED GEANT4 DATA FILES
-----------------------------------------

G4beamline uses Geant4 to perform most of the physics. It inherently needs
data files to drive and control its many physics processes.

To install the data files, run the "g4bldata" program (in the bin directory
under the install directory). It is a GUI program that lets you select which
datasets to download (depends on which physics lists you intend to use).
Note that by default it puts these data files into $HOME/Geant4Data (on
Windows C:\Geant4Data). The other programs look for it there, and if it is
not found they will run g4bldata for you.


INITIAL TESTING
---------------

NOTE: The Root libraries issue a benign error message upon the start of g4bl:
	Error: cannot open file "iostream"  (tmpfile):2:
	*** Interpreter error recovered ***
This happens when their dynamic libraries are attached, so it is difficult to
avoid. Just ignore it for now; g4blgui filters these messages out of the output.


GUI: run the g4blgui program. Click the menu item Tools/DoRegressionTests.

Command-Line: run the g4bltest program (in the bin directory under the install
directory).

Both of these run 113 regression tests. They should all pass, though often
some are omitted (depends on configuration).


RUNNING THE PROGRAM -- GUI
--------------------------

Simply double-click the G4beamline icon. On Windows it is placed on your desktop
and in the Start/G4beamline menu. On Mac OS X it is placed where you dragged 
(copied) it when you installed it. On Linux the bin/g4bl-icon script will create
one on your desktop.

To run the GUI program via the command line:
	source G4beamline-VERSION/bin/g4bl-setup.sh
	# optionally cd where your input file is located
	g4blgui  [input.file]
This requires the X-Windows display on Linux, nothing special on Mac and
Windows. Its opening screen decribes how to use it.

To run the examples, simply push the Browse button, or the File/Open menu item,
navigate to the install directory / examples, and select Example1.g4bl (or 
other *.g4bl file).  Then select the desired viewer (if any), and push the 
Run button.  On Windows, a copy of "G4beamline Examples" is put into
"Documents".


RUNNING THE PROGRAM -- COMMAND-LINE
-----------------------------------

For command-line use (Linux, Mac OS X, and Windows with Cygwin), you must
add G4beamline's bin directory to your PATH:
	source G4beamline-VERSION/bin/g4bl-setup.sh
You can put this into $HOME/.bash_profile. There is also g4bl-setup.csh.

Then cd to whatever directory you plan to use for developing your
simulation(s), and execute:
	g4bl -

NOTE: Root emits a bunch of benign error messages before the "G4BL_DIR=..."
line. Ignore them. The Root team might issue a new version or workaround that
omits them.

After a few seconds it should type (with obvious variations):
    G4BL_DIR=/Users/g4bl/G4beamline-3.02
    G4ABLADATA=/Users/g4bl/Geant4Data/G4ABLA3.0
    G4LEDATA=/Users/g4bl/Geant4Data/G4EMLOW6.41
    G4ENSDFSTATEDATA=/Users/g4bl/Geant4Data/G4ENSDFSTATE1.0
    G4NEUTRONHPDATA=/Users/g4bl/Geant4Data/G4NDL4.5
    G4NEUTRONXSDATA=/Users/g4bl/Geant4Data/G4NEUTRONXS1.4
    G4PIIDATA=/Users/g4bl/Geant4Data/G4PII1.3
    G4SAIDXSDATA=/Users/g4bl/Geant4Data/G4SAIDDATA1.1
    G4LEVELGAMMADATA=/Users/g4bl/Geant4Data/PhotonEvaporation3.1
    G4RADIOACTIVEDATA=/Users/g4bl/Geant4Data/RadioactiveDecay4.2
    G4REALSURFACEDATA=/Users/g4bl/Geant4Data/RealSurface1.0
    G4beamline Process ID 93685
    
    *************************************************************
     g4beamline version: 3.02                        (Jan 26 2016)
                          Copyright : Tom Roberts, Muons, Inc.
                            License : Gnu Public License
                                WWW : http://g4beamline.muonsinc.com
    *************************************************************
     Geant4 version Name: geant4-10-02    (4-December-2015)
                          Copyright : Geant4 Collaboration
                          Reference : NIM A 506 (2003), 250-303
                                WWW : http://cern.ch/geant4
    *************************************************************
    
    geometry                   nPoints=100 printGeometry=0 visual=0
                               tolerance=0.002
    cmd: 

The program is now ready for input. Type "help" to get a short list of
the input-file commands, "help beam" to get help on the beam comand, or 
"help *" for a detailed description of all commands. This is a useful way
to get help on commands when editing your input file(s). Type ^C to exit
back to a shell prompt.


EXAMPLE 1
---------

The first example input file is Example1.g4bl. It is a simple file to track
muons through 1-meter drift spaces into 4 detectors. To visualize its
geometry using the Qt viewer, execute:
	cd G4beamline-VERSION/examples
	g4bl example1.g4bl viewer=best

To run the beam through the geometry, execute:
	g4bl example1.g4bl 
This will generate a Root file named g4beamline.root. To view it do:
	historoot g4beamline.root


OTHER EXAMPLES
--------------

There are other examples in the examples directory.
