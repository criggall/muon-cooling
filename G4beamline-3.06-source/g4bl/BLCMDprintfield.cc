//	BLCMDprintfield.cc
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

#include <stdio.h>
#include <stdarg.h>
#ifndef WIN32
#include <unistd.h>
#endif

#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLGlobalField.hh"
#include "BLFieldMap.hh"
#include "BLMPI.hh"

#ifdef __CYGWIN__
#include "mysnprintf.hh"
#endif

extern void g4bl_exit(int);

/**	class BLCMDprintfield is a command to print E or B fields
 *
 *	If type=print (the default), prints a 2-d array to stdout.
 *
 *	If type=grid or type=cylinder, BLCMDprintfield will write to a file,
 *	in the format of BLFieldMap. This can be used to combine multiple
 *	overlapping magnets/solenoids into a single map, which will be
 *	more efficient in particle tracking. All nonzero field components
 *	are included.
 **/
class BLCMDprintfield : public BLCommand, public BLCallback {
	G4String type;		// "print", "grid", or"cylinder"
	G4int exit;
	// args for type=print
	G4String field;		// "Bx","By","Bz","Ex","Ey","Ez, Btot, Etot"
	G4String layout;	// 2 letters, row,col, of the set {xyzt}
	G4double x;
	G4double y;
	G4double z;
	G4double t;
	G4double drow;
	G4double dcol;
	G4int nrow;
	G4int ncol;
	// args for type=grid
	G4String file;
	G4String comment;
	G4double X0;
	G4double Y0;
	G4double Z0;
	G4int nX;
	G4int nY;
	G4int nZ;
	G4double dX;
	G4double dY;
	G4double dZ;
	// additional args for type=cylinder:
	G4int nR;
	G4double dR;
public:
	BLCMDprintfield();

	~BLCMDprintfield() { }

	G4String commandName() { return "printfield"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();
	
	/// callback() from BLCallback.
	void callback(int type);
private:
	void do_print();
	void do_points();
	void do_grid();
	void do_cylinder();
};

BLCMDprintfield defaultPrintField;	// registers the command, and holds
					// default values of the arguments.

BLCMDprintfield::BLCMDprintfield()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Prints E or B fields, or writes FieldMap file.");
	setDescription("Prints the value of the electromagnetic field components.\n"
		"For type=print, prints one component of the field\n"
		"in a 2-d table.\n"
		"Any coordinate plane can be printed (XY ... ZT).\n"
		"For type=grid or type=cylinder, writes a file in fieldmap\n"
		"format.\n"
		"Global coordinates are used.\n"
		"Units are Tesla for B and MV/meter for E.\n\n"
		"This command is not placed into the geometry.\n\n"
		"NOTE: This command cannot handle time dependency in the "
		"output BLFieldMap file, but can in the printout.\n\n"
		"Note: if you want to plot field vs position or time, the "
		"'fieldntuple' command is probably better, as it is not "
		"limited to the 2-d paper, is easier to use, and lets "
		"you use existing NTuple plotting tools. If you just "
		"want to test a few points, the 'probefield' command lets "
		"you do that interactively.");

	type = "print";
	exit = 0;
	field = "By";
	layout = "zx";
	x = y = z = t = 0.0;
	drow = dcol = 1.0*cm;
	nrow = ncol = 1;
	X0 = Y0 = Z0 = 0.0;
	nX = nY = nZ = 0;
	dX = dY = dZ = 0.0;
	nR = 0;
	dR = 0.0;
}

int BLCMDprintfield::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() > 0) {
		printError("printfield: Invalid argument\n");
		return -1;
	}

	BLCMDprintfield *p = new BLCMDprintfield(defaultPrintField);
	int retval = p->handleNamedArgs(namedArgs);

	// MPI mode: do it only in rank 1.
	// Non-MPI: isRank1() returns true.
	if(BLMPI::isRank1()) {
		// register the object to be called back to do its print.
		// 1 => after reference particle is tracked.
		BLManager::getObject()->registerCallback(p,1);
	} else {
		printf("printfield: performed only in rank 1\n");
	}

	return retval;
}

void BLCMDprintfield::defineNamedArgs()
{
	argString(type,"type","print, grid, or cylinder.");
	argInt(exit,"exit","Set nonzero to exit after printing field");

	if(argMode == HELP) printf("             Arguments for type=print:\n");
	argString(field,"field","The field to print (Bx,By,Bz,Ex,Ey,Ez,Btot,Etot).");
	argString(layout,"layout","Layout (RowCol) - 2 chars 'AB' each of {xyzt}.");
	argDouble(x,"x","The starting value of x (mm).");
	argDouble(y,"y","The starting value of y (mm).");
	argDouble(z,"z","The starting value of z (mm).");
	argDouble(t,"t","The starting value of time (ns).");
	argDouble(drow,"drow","The incr between points in each row (mm|ns).");
	argDouble(dcol,"dcol","The incr between points in each column (mm|ns).");
	argInt(nrow,"nrow","The number of rows.");
	argInt(ncol,"ncol","The number of columns.");

	if(argMode == HELP) printf("             Arguments for type=grid:\n");
	argString(file,"file","Filename to write fieldmap to.");
	argString(comment,"comment","Comment for fieldmap.");
	argDouble(X0,"X0","Initial value of X (mm, default=0).");
	argDouble(Y0,"Y0","Initial value of Y (mm, default=0).");
	argDouble(Z0,"Z0","Initial value of Z (mm, default=0).");
	argInt(nX,"nX","Number of points in X.");
	argInt(nY,"nY","Number of points in Y.");
	argInt(nZ,"nZ","Number of points in Z.");
	argDouble(dX,"dX","Interval in X between points (mm).");
	argDouble(dY,"dY","Interval in Y between points (mm).");
	argDouble(dZ,"dZ","Interval in Z between points (mm).");

	if(argMode == HELP) printf("             Arguments for type=cylinder:\n");
	argString(file,"file","Filename to write fieldmap to.");
	argString(comment,"comment","Comment for fieldmap.");
	argDouble(Z0,"Z0","Initial value of Z (mm, default=0).");
	argInt(nR,"nR","Number of points in R.");
	argDouble(dR,"dR","Interval in R between points (mm).");
	argInt(nZ,"nZ","Number of points in Z.");
	argDouble(dZ,"dZ","Interval in Z between points (mm).");
}

void BLCMDprintfield::callback(int callback_type)
{
	G4cout.flush();

	// get 1 field value, so warning is printed here, not in loop
	BLGlobalField *gf = BLGlobalField::getObject();
	G4double pos[4], field[6];
	pos[0]=pos[1]=pos[2]=pos[3]=0.0;
	gf->GetFieldValue(pos,field);

	if(type == "print")
		do_print();
	else if(type == "points")
		do_points();
	else if(type == "grid")
		do_grid();
	else if(type == "cylinder")
		do_cylinder();
	else
		printError("printfield: invalid type '%s'\n",type.c_str());

	fflush(stdout);

	if(exit)
		G4Exception("printfield","Exiting",FatalException,"");
}

void BLCMDprintfield::do_print()
{
	int index=-1; 
	if(field == "Bx") index = 0;
	else if(field == "By") index = 1;
	else if(field == "Bz") index = 2;
	else if(field == "Ex") index = 3;
	else if(field == "Ey") index = 4;
	else if(field == "Ez") index = 5;
	else if(field == "Btot") index = 6;
	else if(field == "Etot") index = 7;
	else printError("printfield: invalid field '%s'\n",field.c_str());

	int irow=-1,icol=-1;
	char rowname, colname;
	rowname = layout.c_str()[0];
	colname = layout.c_str()[1];
	switch(rowname) {
	case 'x': case 'X':	irow = 0;	break;
	case 'y': case 'Y':	irow = 1;	break;
	case 'z': case 'Z':	irow = 2;	break;
	case 't': case 'T':	irow = 3;	break;
	}
	switch(colname) {
	case 'x': case 'X':	icol = 0;	break;
	case 'y': case 'Y':	icol = 1;	break;
	case 'z': case 'Z':	icol = 2;	break;
	case 't': case 'T':	icol = 3;	break;
	}

	if(index < 0 || irow < 0 || icol < 0) {
		printError("printfield: Invalid field name or layout\n");
		return;
	}

	BLGlobalField *gf = BLGlobalField::getObject();

	G4double pos[4];
	pos[0]=x, pos[1]=y, pos[2]=z, pos[3]=t;
	G4double col0 = pos[icol];

	printf("\n%s (%s):  ",field.c_str(),(index<=2||index==6 ? "Tesla" : "MV/meter"));
	G4double unit = (index<=2||index==6 ? tesla : megavolt/meter);
	for(int i=0; i<4; ++i) {
		if(i == irow || i == icol) continue;
		printf("%c=%.3f ","xyzt"[i],pos[i]);
	}
	printf("\n           ");
	for(int col=0; col<ncol; ++col) {
		char tmp[16];
		sprintf(tmp,"%c=%.3f",colname,pos[icol]);
		printf("%11.11s",tmp);
		pos[icol] += dcol;
	}
	printf("\n");
	for(int row=0; row<nrow; ++row) {
		printf("%c=%-8.3f ",rowname,pos[irow]);
		pos[icol] = col0;
		for(int col=0; col<ncol; ++col) {
			G4double field[6];
			gf->GetFieldValue(pos,field);
			if(index == 6)
				printf(" %10.4f",sqrt(field[0]*field[0]+field[1]*field[1]+field[2]*field[2])/unit);
			else if(index == 7)
				printf(" %10.4f",sqrt(field[3]*field[3]+field[4]*field[4]+field[5]*field[5])/unit);
			else
				printf(" %10.4f",field[index]/unit);
			pos[icol] += dcol;
		}
		printf("\n");
		pos[irow] += drow;
	}

	printf("\n");

}

void BLCMDprintfield::do_points()
{
	int irow=-1,icol=-1;
	char rowname, colname;
	rowname = layout.c_str()[0];
	colname = layout.c_str()[1];
	switch(rowname) {
	case 'x': case 'X':	irow = 0;	break;
	case 'y': case 'Y':	irow = 1;	break;
	case 'z': case 'Z':	irow = 2;	break;
	case 't': case 'T':	irow = 3;	break;
	}
	switch(colname) {
	case 'x': case 'X':	icol = 0;	break;
	case 'y': case 'Y':	icol = 1;	break;
	case 'z': case 'Z':	icol = 2;	break;
	case 't': case 'T':	icol = 3;	break;
	}

	if(irow < 0 || icol < 0) {
		printError("printfield: Invalid field name or layout\n");
		return;
	}

	BLGlobalField *gf = BLGlobalField::getObject();

	G4double pos[4];
	pos[0]=x, pos[1]=y, pos[2]=z, pos[3]=t;
	G4double col0 = pos[icol];

	printf("# B field: x y z Bx By Bz\n");
	for(int row=0; row<nrow; ++row) {
		pos[icol] = col0;
		for(int col=0; col<ncol; ++col) {
			G4double field[6];
			gf->GetFieldValue(pos,field);
			printf("%.1f %.1f %.1f %.4f %.4f %.4f\n",
				pos[0],pos[1],pos[2],
				field[0]/tesla,field[1]/tesla,field[2]/tesla);
			pos[icol] += dcol;
		}
		pos[irow] += drow;
	}

	if(exit) {
		printf("printfield: exit\n");
		g4bl_exit(0);
	}
}

void BLCMDprintfield::do_grid()
{
	if(file == "" || nX <= 0 || nY <= 0 || nZ <= 0 || 
	   dX <= 0.0 || dY <= 0.0 || dZ <= 0.0) {
		printError("printfield type=grid: invalid arguments\n");
		return;
	}
	BLFieldMap fm;
	fm.createGridMap(X0,Y0,Z0,dX,dY,dZ,nX,nY,nZ,BLGlobalField::getObject());
	fm.writeFile(file, comment);
}

void BLCMDprintfield::do_cylinder()
{
	if(file == "" || nR <= 0 || nZ <= 0 || dR <= 0.0) {
		printError("printfield type=cylinder: invalid arguments\n");
		return;
	}
	BLFieldMap fm;
	fm.createCylinderMap(Z0,dR,dZ,nR,nZ, BLGlobalField::getObject());
	fm.writeFile(file, comment);
}

