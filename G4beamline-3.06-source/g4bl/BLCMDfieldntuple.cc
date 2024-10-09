//	BLCMDfieldntuple.cc
/*
This source file is part of G4beamline, http://g4beamline.muonsinc.com
Copyright (C) 2002-2013 by Tom Roberts, all rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
*/

#include <vector>

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLGlobalField.hh"
#include "BLNTuple.hh"
#include "BLMPI.hh"

/**	class BLCMDfieldntuple is a command to generate an NTuple from field
 *	values at specified points.
 *
 **/
class BLCMDfieldntuple : public BLCommand, public BLCallback {
	G4String name;
	G4String category;
	G4String format;
	G4String filename;
	G4String x;
	G4String y;
	G4String z;
	G4String t;
	G4int exit;
	class Loop {
		double &coord;
		std::vector<double> value;
		unsigned index;
	public:
		Loop(double &_coord, G4String s);
		void start() { index=0; }
		bool next() {
			if(index >= value.size()) return false;
			coord = value[index++];
			return true;
		}
	};
public:
	BLCMDfieldntuple();

	~BLCMDfieldntuple() { }

	BLCMDfieldntuple(const BLCMDfieldntuple& r);

	G4String commandName() { return "fieldntuple"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();
	
	/// callback() from BLCallback.
	void callback(int type);
};

BLCMDfieldntuple defaultFieldntuple;	// registers the command, and holds
					// default values of the arguments.
BLCMDfieldntuple::BLCMDfieldntuple() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Generates an NTuple from B and E fields at specified points.");
	setDescription("Intended primarily for debugging. "
		"This command makes it easy to plot fields as a function of "
		"position and time, using existing NTuple plotting tools. "
		"Outputs x,y,z,t,Bx,By,Bz,Ex,Ey,Ez into an NTuple. "
		"Units are mm, ns, Tesla, and MegaVolts/meter. "
		"Runs after the reference particle is tracked. "
		"Only global coordinates are used.\n\n"
		"The single positional argument is the name of the NTuple. "
		"Named arguments {x,y,z,t} are of two forms specifying "
		"coordinate values: "
		"x=Xmin,Xmax,dX or x=X1:X2:X3:... which generate the obvious "
		"loops (single value is OK). Expressions can be used. "
		"Omitted coordinates are held fixed at 0.0.\n\n"
		"This command is not placed into the geometry.");
	name = "";
	category = "";
	format = "";
	filename = "";
	x = "0.0";
	y = "0.0";
	z = "0.0";
	t = "0.0";
	exit = 0;
}

BLCMDfieldntuple::BLCMDfieldntuple(const BLCMDfieldntuple& r) : BLCommand(), BLCallback()
{
	name = r.name;
	category = r.category;
	format = r.format;
	filename = r.filename;
	x = r.x;
	y = r.y;
	z = r.z;
	t = r.t;
	exit = r.exit;
}

int BLCMDfieldntuple::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("fieldntuple: Invalid name\n");
		return -1;
	}

	BLCMDfieldntuple *p = new BLCMDfieldntuple(defaultFieldntuple);
	p->name = argv[0];
	int retval = p->handleNamedArgs(namedArgs);

	// MPI mode: generate the NTuple only in rank 1;
	// non-MPI: isRank1() returns true.
	if(BLMPI::isRank1()) {
		// register the object to be called back to do its work.
		// 1 => after reference particle is tracked.
		BLManager::getObject()->registerCallback(p,1);
	} else {
		printf("fieldntuple: only performed in rank 1\n");
	}

	return retval;
}

void BLCMDfieldntuple::defineNamedArgs()
{
	argString(category,"category","The category of the NTuple.");
	argString(format,"format","NTuple format (see above for list).");
	argString(filename,"filename","Name of file.");
	argString(x,"x","Loop for x values: Xmin,Xmaz,dX or X1:X2:... (mm)");
	argString(y,"y","Loop for y values: Ymin,Ymaz,dY or Y1:Y2:... (mm)");
	argString(z,"z","Loop for z values: Zmin,Zmaz,dZ or Z1:Z2:... (mm)");
	argString(t,"t","Loop for t values: Tmin,Tmaz,dT or T1:T2:... (ns)");
	argInt(exit,"exit","Set nonzero to exit after generating NTuple (0).");
	argString(filename,"file","Synonym for filename.");
}

void BLCMDfieldntuple::callback(int callback_type)
{
	BLGlobalField *gf = BLGlobalField::getObject();

	double point[4], field[6];
	Loop loopX(point[0],x);
	Loop loopY(point[1],y);
	Loop loopZ(point[2],z);
	Loop loopT(point[3],t);

	BLNTuple *ntuple=BLNTuple::create(format,category,name,
					"x:y:z:t:Bx:By:Bz:Ex:Ey:Ez",filename);

	for(loopX.start(); loopX.next(); ) {
	    for(loopY.start(); loopY.next(); ) {
		for(loopZ.start(); loopZ.next(); ) {
		    for(loopT.start(); loopT.next(); ) {
			gf->GetFieldValue(point,field);
			double data[10];
			data[0] = point[0]/mm;
			data[1] = point[1]/mm;
			data[2] = point[2]/mm;
			data[3] = point[3]/ns;
			data[4] = field[0]/tesla;
			data[5] = field[1]/tesla;
			data[6] = field[2]/tesla;
			data[7] = field[3]/(megavolt/meter);
			data[8] = field[4]/(megavolt/meter);
			data[9] = field[5]/(megavolt/meter);
			ntuple->appendRow(data,10);
		    }
		}
	    }
	}

	ntuple->flush();

	if(exit != 0)
		G4Exception("fieldntuple","Exiting",FatalException,"");
}

BLCMDfieldntuple::Loop::Loop(double &_coord, G4String s) : coord(_coord)
{
	value = getList(s,",");
	if(value.size() == 3) {
		double a = value[0];
		double b = value[1];
		double c = value[2];
		value.clear();
		if(c > 0.0 && b >= a) {
			for(double v=a; v<=b; v+= c)
				value.push_back(v);
		}
	} else {
		value = getList(s,":");
	}
	if(value.size() == 0) {
		printError("Invalid coordinate loop.");
		value.push_back(0.0);
	}
	index = 0;
	coord = value[0];
}

