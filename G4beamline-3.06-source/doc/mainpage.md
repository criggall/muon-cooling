\mainpage g4beamline

g4beamline is a program to run geant4 simulations of beamlines. It is
designed to be easy to use, extremely flexible, and easily extensible.

<h2>Extending the g4beamline program</h2>

Additional elements and commands can easily be added to the g4beamline program
by simply deriving additional classes from BLElement or BLCommand. The
source has lots of examples of each. The basic idea is that for each
command or element the .cc file defines a global default instance, and
the default constructor registers the command with BLCommand so it can
be used in the input file. Virtual functions from BLCommand define the
parameters of the command (see BLCommand for details). New elements can
be placed (via the "place" command) as usual.

Because of the rich set of callbacks (see below), a new command can do just
about anything that is possible in the program. Commands generally create
a new instance of the class, and register it for the appropriate callbacks
to implement the functionality required. This can include Geant4 callbacks
and functions such as registering new particles and physics processes.

Most new features in G4beamline are implemented as new commands. Occasionally
a new feature will require modifications to the infrastructure classes, or 
even new infrastructure classes, but that is rather rare. Multiple new
commands can be developed and tested in parallel by different people in their
private source trees; once one is ready, its file is placed into the official
source tree and it is automatically included in the next official build.
This can be done aggressively, because if the new command is not used in the
input file, its code will not be executed (except for registering its command).
So even if the new command has buggy code, it won't affect users that don't 
use it.

<h2>Implementation Notes</h2>

The commands for an input file are all implemented as individual classes in
individual .cc files, derived from BLCommand. The class BLCommand
handles reading the input file, parsing commands, parsing
command arguments (both positional and named), and implementing first-character
commands (/, *, #, and !). Derived commands simply declare their arguments (via 
calls to argString(), argInt(), and argDouble()), and handleNamedArgs()
manages putting the values from the command line into the class variables.

Conventionally, the class name, filename, and command name are interrelated:
<pre>
	command name:	example
	class name:	BLCMDexample
	file name:	BLCMDexample.cc
</pre>
Commands do not require a .hh file, as the class is used only in the
command's single .cc file. Classes used in more than one command should
be implemented as infrastructure classes.

All .cc files in the <i>source/src</i> directory are included in the build,
so their global static constructors will be executed before <i>main()</i>
is entered. As those constructors register the default instance of each
command with BLCommand, no list of commands is required -- simply placing
the source file into that directory is sufficient.

Beamline elements are also implemented as individual classes in individual
.cc files, derived from BLElement (which is derived from BLCommand, so
every element includes a command to create an instance of the element).
The place command places instances of a given element into the current group.

The class BLManager interfaces to BLRunManager, and implements <i>all</i>
of the user customization classes (G4UserSteppingAction, 
G4UserDetectorConstruction, etc.). Element classes can register with BLManager
to have a UserSteppingAction() function called for every reference-particle
or beam step referring to a specific G4VPhysicalVolume. For a simple example
of this see BLCMDvirtualdetector.cc.

<h2>Callback Summary</h2>
The BLManager provides callback services to command and element classes; the
latter register with the former, specifying when to be called back.
Some of these are extensions of Geant4 callbacks, some are new in G4beamline.
Note the Geant4 callbacks occur in every state, so the code should check the 
BLManager state.

For a simulation run (i.e. viewer=none), the order of callbacks is:
<table border=1>
<tr><th>Callback</th><th>Description</th></tr>
<tr>
<td>Reading Input File</td><td>Commands execute, elements construct their objects. For any command class, reading its command generates a call to <i>command()</i>, while for elements the place command generates a call to <i>construct()</i>. Most element commands create a new class instance in <i>command()</i> which holds the values of its arguments; the heavy work is performed in <i> construct()</i>, such as creating solids, logical and physical volumes, and placing them into the Geant4 world.</td>
</tr>
<tr>
<td>callback(0)</td><td>Callback pre-Tune particle.</td>
</tr>
<tr>
<td>-</td><td>Tune particle is tracked, including BeginOfRunAction, BeginOfEventAction, PreUserTrackingAction (multiple times), UserSteppingAction (multiple times), ZSteppingAction (multiple times), StackingAction (multiple times), PostUserTrackingAction, EndOfEventAction, EndOfRunAction.</td>
</tr>
<tr>
<td>-</td><td>Reference particle is tracked, including BeginOfRunAction, BeginOfEventAction, PreUserTrackingAction, UserSteppingAction (multiple times), ZSteppingAction (multiple times), StackingAction, PostUserTrackingAction, EndOfEventAction, EndOfRunAction.</td>
</tr>
<tr>
<td>callback(1)</td><td>Callback post-Reference particle.</td>
</td>
<tr>
<td>BeginOfRunAction</td><td>Start of a run tracking beam.</td>
</tr>
<tr>
<td>BeginOfEventAction</td><td>Just before each event is processed.</td>
</tr>
<tr>
<td>PreUserTrackingAction</td><td>Just before each track is tracked.</td>
</tr>
<tr>
<td>UserSteppingAction</td><td>During tracking, each 'physics' step. This is the most common callback, and can be conditioned on the step being in a specific physical volume, and on the state (Tune, Reference, Beam).</td>
</tr>
<tr>
<td>ZSteppingAction</td><td>During tracking, each 'physics' step at a specific Z value (centerline coordinates).</td>
</tr>
<tr>
<td>StackingAction</td><td>During tracking, when placing secondary track onto stack.</td>
<tr>
<td>PostUserTrackingAction</td><td>Just after each track is tracked.</td>
</tr>
<tr>
<td>EndOfEventAction</td><td>Just after each event is processed.</td>
</tr>
<tr>
<td>EndOfRunAction</td><td>End of a run tracking beam.</td>
</tr>
</tr>
<tr>
<td>callback(2)</td><td>Callback post-beam tracking.</td>
</td>
</table>

Note that the Tune particle consists of a run with one event containing one or more tracks (depends on how often the Tune particle must be re-tracked during tuning). The Reference particle consists of a run with one event containing one track. So the corresponding Geant4 callbacks happen for these two runs as well as for the run containing all of the beam tracks.

For visualization, callback(4) is called after callback(1) and after the selected viewer is set up. For each image, a run is simulated with visualization of its tracks enabled.

For special situations, such as collective tracking, callback(3) is called just after callback(1), and the program exits after all registered callback(3)-s return.

There are additional callbacks that are not related to runs, events, or tracking:
<table border=0>
<tr> <td>BLCommand</td><td>Registers command names and maps them to implementation classes.</td> </tr>
<tr><td>BLCollectiveComputation</td><td>Registered with BLRunManager to perform collective computations.</td></tr>
<tr><td>BLManager</td><td>In addition to the above callbacks, it also registers the BLPhysics object, BLBeam objects, and BLUserCode objects.</td></tr>
<tr><td>BLNTuple</td><td>Registers NTuple Handlers to implement the different types of NTuples, and also registers callbacks to be called by the <i>appendRow()</i> function of specific NTuples.</td></tr>
<tr><td>BLGlobalField</td><td>Registers individual instances of BLElementField that implement the electromagnetic field of an individual element. Handles field overlaps, and uses bounding boxes for efficiency.</td></tr>
</table>


