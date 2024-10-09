//	BLFieldMap.cc
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


#include <fstream>
#include <string.h>
#include <ctype.h>

#include "G4ElectroMagneticField.hh"

#include "BLAssert.hh"
#include "BLFieldMap.hh"
#include "BLCommand.hh"		// command-parsing routines only
#include "BTSpline1D.hh"

/**	class InputFile handles input file I/O 
 **/
class InputFile {
	std::ifstream in;
	G4String name;
	char *line;
	int maxline;
	int lineno;
	bool repeatline;
public:
	InputFile(G4String filename, G4int _maxline=1024) : in(), name() {
		name = filename;
		in.open(name.c_str());
		maxline = _maxline;
		line = new char[maxline];
		BLAssert(line != 0);
		lineno = 0;
		repeatline = false;
	}
	~InputFile() { close(); }
	const char *filename() { return name.c_str(); }
	bool good() { return in.good(); }
	int linenumber() { return lineno; }
	char *getline() {		// skips blank lines and comments
		if(repeatline) {
			repeatline = false;
			return line;
		}
		while(in.good()) {
			line[0] = '\0';
			in.getline(line,maxline);
			++lineno;
			if(line[0] == '*')
				printf("%s\n",line);
			if(line[0] == '\0' || line[0] == '#' || line[0] == '*')
				continue;
			return line;
		}
		return 0;
	}
	void setMaxline(int _maxline) {
		maxline = _maxline;
		BLAssert(maxline >= 80);
		if(line) delete line;
		line = new char[maxline];
		BLAssert(line != 0);
	}
	void close() { if(line) { in.close(); delete line; line = 0; } }
	void repeatLine() { repeatline = true; }
};

/**	class FieldMapImpl -- base class for a FieldMap implementation
 **/
class FieldMapImpl {
public:
	FieldMapImpl() { }
	virtual ~FieldMapImpl() { }
	virtual void getFieldValue(const G4double local[4], G4double field[6])
								const = 0;
	virtual bool handleCommand(InputFile &in, BLArgumentVector &argv,
			BLArgumentMap &namedArgs) = 0;
	virtual void getBoundingPoint(int i, G4double point[4]) = 0;
	virtual bool hasB() = 0;
	virtual bool hasE() = 0;
	bool readBlock(InputFile &in, float *values, int nRows, int nCols,
							G4double units);
	virtual bool writeFile(FILE *f) = 0;
};

/**	class GridImpl -- class for a Grid FieldMap implementation
 **/
class GridImpl : public FieldMapImpl {
	G4int nX;
	G4int nY;
	G4int nZ;
	G4double dX;
	G4double dY;
	G4double dZ;
	G4double X0;
	G4double Y0;
	G4double Z0;
	G4double tolerance;
	float *mapBx;
	float *mapBy;
	float *mapBz;
	float *mapEx;
	float *mapEy;
	float *mapEz;
	bool extendX;
	bool extendY;
	bool extendZ;
	int extendXbits;
	int extendYbits;
	int extendZbits;
public:
	GridImpl(BLArgumentVector &argv, BLArgumentMap &namedArgs);
	~GridImpl();
	void getFieldValue(const G4double local[4], G4double field[6]) const;
	bool handleCommand(InputFile &in, BLArgumentVector &argv,
			BLArgumentMap &namedArgs);
	virtual void getBoundingPoint(int i, G4double point[4]);
	virtual bool hasB() { return mapBx!=0 || mapBy!=0 || mapBz!=0; }
	virtual bool hasE() { return mapEx!=0 || mapEy!=0 || mapEz!=0; }
	int bits(G4String s);
	const char *bits2str(int v);
	bool setField(G4double X, G4double Y, G4double Z, G4double Bx,
		G4double By, G4double Bz, G4double Ex, G4double Ey,
		G4double Ez, int linenumber);
	virtual bool writeFile(FILE *f);
};

/**	class CylinderImpl -- class for a Cylinder FieldMap implementation
 **/
class CylinderImpl : public FieldMapImpl {
	G4int nR;
	G4int nZ;
	G4double dR;
	G4double dZ;
	G4double Z0;
	G4double tolerance;
	float *mapBz;
	float *mapBr;
	float *mapBphi;
	float *mapEr;
	float *mapEz;
	bool extendZ;
	float extendBrFactor, extendBzFactor;
	float extendErFactor, extendEzFactor;
public:
	CylinderImpl(BLArgumentVector &argv, BLArgumentMap &namedArgs);
	~CylinderImpl();
	void getFieldValue(const G4double local[4], G4double field[6]) const;
	bool handleCommand(InputFile &in, BLArgumentVector &argv,
			BLArgumentMap &namedArgs);
	virtual void getBoundingPoint(int i, G4double point[4]);
	virtual bool hasB() { return mapBz != 0 || mapBr != 0 || mapBphi != 0; }
	virtual bool hasE() { return mapEz != 0 || mapEr != 0; }
	bool setField(G4double R, G4double Z, G4double Bz, G4double Br,
		G4double Bphi, G4double Er, G4double Ez, int linenumber);
	virtual bool writeFile(FILE *f);
};

/**	class TimeImpl -- class implementation of a time dependence.
 */
class TimeImpl {
	BTSpline1D B;
	BTSpline1D E;
	G4double tMin;
	G4double tMax;
	G4double period;
public:
	/// Constructor.
	/// n is # points in the arrays, t[] is the list of time values (ns),
	/// b[] and e[] are the factors for B field and E field.
	/// t[], b[], and e[] are copied so the caller can
	/// delete them after calling this routine.
	TimeImpl(int n, G4double t[], G4double b[], G4double e[]=0) :
						B(n,t,b), E(n,t,(e?e:b)) {
		tMin = t[0];
		tMax = t[n-1];
		period = -1.0;
	}
	void setPeriod(G4double v) { period = v; }
	G4double factorB(G4double t) { 
		if(period > 0.0) {
			t -= floor(t/period)*period;
		} else {
			if(t < tMin) t = tMin;
			if(t > tMax) t = tMax;
		}
		return B(t);
	}
	G4double factorE(G4double t) { 
		if(period > 0.0) {
			t -= floor(t/period)*period;
		} else {
			if(t < tMin) t = tMin;
			if(t > tMax) t = tMax;
		}
		return E(t);
	}
	static TimeImpl *readTime(InputFile &in, BLArgumentVector argv,
						BLArgumentMap namedArgs);
};

/// argDouble handles an argument with a value of G4double
static void argDouble(G4double &var, const char *name, BLArgumentMap &namedArgs)
{
	if(namedArgs.count(name) == 0) return;
	char *q;
	const char *p=namedArgs[name].c_str();
	G4double v = strtod(p,&q);
	if(p == q || *q != '\0') {
		BLCommand::printError("Invalid value for '%s' in BLFieldMap",
								name);
		return;
	}
	var = v;
}


/// argInt handles an argument with a value of G4int
static void argInt(G4int &var, const char *name, BLArgumentMap &namedArgs)
{
	if(namedArgs.count(name) == 0) return;
	char *q;
	const char *p=namedArgs[name].c_str();
	G4int v = strtol(p,&q,0);
	if(p == q || *q != '\0') {
		BLCommand::printError("Invalid value for '%s' in BLFieldMap",
								name);
		return;
	}
	var = v;
}

/// d2string converts a double to a G4String
static G4String d2string(G4double v)
{
	char tmp[32];
	sprintf(tmp,"%.8g",v);
	return G4String(tmp);
}

/// i2string converts an int to a G4String
static G4String i2string(int v)
{
	char tmp[32];
	sprintf(tmp,"%d",v);
	return G4String(tmp);
}

BLFieldMap::BLFieldMap()
{
	maxline = 1024;
	current = 1.0;
	gradient = 1.0;
	normB = 1.0;
	normE = 1.0;
	impl = 0;
	time = 0;
}

BLFieldMap::~BLFieldMap()
{
	if(impl) delete impl;
	impl = 0;
	if(time) delete time;
	time = 0;
}

void BLFieldMap::getFieldValue(const G4double local[4], G4double field[6],
                        G4double _current, G4double _gradient)
{
	if(!impl)
		throw "BLFieldMap::getFieldValue called, no implementation";

	G4double thisField[6];
	impl->getFieldValue(local,thisField);
	G4double timeB=1.0, timeE=1.0;
	if(time) {
		timeB = time->factorB(local[3]);
		timeE = time->factorE(local[3]);
	}
	field[0] = thisField[0] * normB * timeB * _current/current;
	field[1] = thisField[1] * normB * timeB * _current/current;
	field[2] = thisField[2] * normB * timeB * _current/current;
	field[3] = thisField[3] * normE * timeE * _gradient/gradient;
	field[4] = thisField[4] * normE * timeE * _gradient/gradient;
	field[5] = thisField[5] * normE * timeE * _gradient/gradient;
}

//KBB-8mar11: - avoid rescaling & time dependence (used in BLCMDpillbox)
//              current,gradient just a scaling (default=1)
void BLFieldMap::getFieldValueNoTimeNoScale(const G4double local[4], G4double field[6])
{
	if(!impl)
		throw "BLFieldMap::getFieldValueNoTimeNoScale called, no implementation";

	impl->getFieldValue(local,field);
}

bool BLFieldMap::readFile(G4String filename)
{
	InputFile in(filename,maxline);
	if(!in.good()) {
		BLCommand::printError("BLFieldMap Cannot open file '%s'",
				in.filename());
		return false;
	}
	printf("BLFieldMap: reading file '%s'\n",in.filename());

	bool retval = true;

	char *line;
	while((line=in.getline()) != 0) {
		BLArgumentVector argv;
		BLArgumentMap namedArgs;
		if(BLCommand::parseArgs(line,argv,namedArgs) < 0)
			goto invalid;
		if(argv[0] == "") continue;
		if(argv[0] == "param") {
			argInt(maxline,"maxline",namedArgs);
			argDouble(current,"current",namedArgs);
			argDouble(gradient,"gradient",namedArgs);
			argDouble(normB,"normB",namedArgs);
			argDouble(normE,"normE",namedArgs);
			in.setMaxline(maxline);
		} else if(argv[0] == "grid") {
			if(impl) goto invalid;
			impl = new GridImpl(argv,namedArgs);
		} else if(argv[0] == "cylinder") {
			if(impl) goto invalid;
			impl = new CylinderImpl(argv,namedArgs);
		} else if(argv[0] == "time") {
			if(time) goto invalid;
			time = TimeImpl::readTime(in,argv,namedArgs);
			if(!time) goto invalid;
		} else if(impl) {
			if(!impl->handleCommand(in,argv,namedArgs))
				goto invalid;
		} else {
invalid:		BLCommand::printError("BLFieldMap file '%s' line %d invalid command '%s'\n",
					in.filename(),in.linenumber(),
					argv[0].c_str());
			retval = false;
			break;
		}
	}
	in.close();

	return retval;
}


void BLFieldMap::getBoundingPoint(int i, G4double point[4])
{
	impl->getBoundingPoint(i,point);
}

bool BLFieldMap::hasB()
{
	return impl->hasB();
}

bool BLFieldMap::hasE()
{
	return impl->hasE();
}


TimeImpl *TimeImpl::readTime(InputFile &in, BLArgumentVector argv,
						BLArgumentMap namedArgs)
{
	G4double period=-1.0;
	argDouble(period,"period",namedArgs);

	std::vector<G4double> T;
	std::vector<G4double> B;
	std::vector<G4double> E;

	char *line;
	while((line=in.getline()) != 0) {
		if(isalpha(line[0])) {
			in.repeatLine();
			break;
		}
		G4double t,b,e;
		int j = sscanf(line,"%lf%lf%lf",&t,&b,&e);
		if(j < 2 || j > 3) return 0;
		if(j == 2) e = b;
		T.push_back(t);
		B.push_back(b);
		E.push_back(e);
	}

	if(T.size() > 0) {
		TimeImpl *impl =
			new TimeImpl(T.size(),&T[0],&B[0],&E[0]);
		impl->setPeriod(period);
		return impl;
	}

	return 0;
}

bool FieldMapImpl::readBlock(InputFile &in, float *values, int nRows, int nCols,
						G4double units)
{
	while(nRows-- > 0) {
		char *line = in.getline();
		if(!line) return false;
		char *p=line;
		for(int i=0; i<nCols; ++i) {
			while(isspace(*p)) ++p;
			if(*p == '\0') return false;
			*values++ = strtod(p,&p) * units;
			if(*p == ',') ++p;
		}
	}
	return true;
}

bool BLFieldMap::createGridMap(G4double X0, G4double Y0, G4double Z0, 
	G4double dX, G4double dY, G4double dZ, int nX, int nY, int nZ,
	G4ElectroMagneticField *emField)
{
	if(impl) {
		delete impl;
		impl = 0;
	}
	maxline = 128;
	current = 1.0;
	gradient = 1.0;
	normB = 1.0;
	normE = 1.0;

	BLArgumentVector argv;
	BLArgumentMap args;
	args["X0"] = d2string(X0);
	args["Y0"] = d2string(Y0);
	args["Z0"] = d2string(Z0);
	args["dX"] = d2string(dX);
	args["dY"] = d2string(dY);
	args["dZ"] = d2string(dZ);
	args["nX"] = i2string(nX);
	args["nY"] = i2string(nY);
	args["nZ"] = i2string(nZ);

	GridImpl *grid = new GridImpl(argv,args);
	impl = grid;

	bool retval = true;
	G4double pos[4], field[6];
	pos[3] = 0.0;
	for(int i=0; i<nX; ++i) {
		pos[0] = X0 + i*dX;
		for(int j=0; j<nY; ++j) {
			pos[1] = Y0 + j*dY;
			for(int k=0; k<nZ; ++k) {
				pos[2] = Z0 + k*dZ;
				emField->GetFieldValue(pos,field);
				if(!grid->setField(pos[0],pos[1],pos[2],
					    field[0],field[1],field[2],
					    field[3],field[4],field[5],0))
					retval = false;
			}
		}
	}

	return retval;
}

bool BLFieldMap::createCylinderMap(G4double Z0, G4double dR,  G4double dZ,
	int nR, int nZ, G4ElectroMagneticField *emField)
{
	if(impl) {
		delete impl;
		impl = 0;
	}
	maxline = 128;
	current = 1.0;
	gradient = 1.0;
	normB = 1.0;
	normE = 1.0;

	BLArgumentVector argv;
	BLArgumentMap args;
	args["Z0"] = d2string(Z0);
	args["dR"] = d2string(dR);
	args["dZ"] = d2string(dZ);
	args["nR"] = i2string(nR);
	args["nZ"] = i2string(nZ);

	CylinderImpl *cyl = new CylinderImpl(argv,args);
	impl = cyl;

	// use the Y=0 plane for the R,Z plane
	bool retval = true;
	G4double pos[4], field[6];
	pos[3] = 0.0;
	for(int i=0; i<nR; ++i) {
		pos[0] = 0.0 + i*dR;
		pos[1] = 0.0;
		for(int k=0; k<nZ; ++k) {
			pos[2] = Z0 + k*dZ;
			emField->GetFieldValue(pos,field);
			if(!cyl->setField(pos[0],pos[2],field[2],field[0],
					    field[1],field[3],field[5],0))
				retval = false;
		}
	}

	return false;
}

bool BLFieldMap::writeFile(G4String filename, G4String comment)
{
	if(!impl) return false;

	FILE *f = fopen(filename.c_str(),"r");
	if(f) {
		fclose(f);
		G4Exception("BLFieldMap","Output File Exists",FatalException,
							filename.c_str());
	}

	f = fopen(filename,"w");
	if(!f) {
		fprintf(stderr,"BLFieldMap::writeFile CANNOT WRITE file '%s'\n",
				filename.c_str());
		return false;
	}

	fprintf(f,"# %s\n",comment.c_str());
	fprintf(f,"param maxline=%d current=%g gradient=%g normB=%g normE=%g\n",
		maxline,current,gradient,normB,normE);

	bool retval = impl->writeFile(f);

	fclose(f);
	return retval;
}

bool BLFieldMap::createTimeDependence(int n, G4double t[], G4double b[],
					G4double e[], G4double period)
{
	if(time) return false;
	time = new TimeImpl(n,t,b,e);
	if(period > 0.0)
		time->setPeriod(period);
	return true;
}

bool BLFieldMap::getTimeFactor(G4double t, G4double *b, G4double *e)
{
	if(b) *b = (time ? time->factorB(t) : 1.0);
	if(e) *e = (time ? time->factorE(t) : 1.0);
	return true;
}







GridImpl::GridImpl(BLArgumentVector &argv, BLArgumentMap &namedArgs)
							: FieldMapImpl()
{
	nX = 2;
	nY = 2;
	nZ = 2;
	dX = 10.0*mm;
	dY = 10.0*mm;
	dZ = 10.0*mm;
	X0 = 0.0;
	Y0 = 0.0;
	Z0 = 0.0;
	tolerance = 0.01*mm;
	mapBx = 0;
	mapBy = 0;
	mapBz = 0;
	mapEx = 0;
	mapEy = 0;
	mapEz = 0;
	extendX = false;
	extendY = false;
	extendZ = false;
	extendXbits = 0;
	extendYbits = 0;
	extendZbits = 0;
	argInt(nX,"nX",namedArgs);
	argInt(nY,"nY",namedArgs);
	argInt(nZ,"nZ",namedArgs);
	argDouble(dX,"dX",namedArgs);
	argDouble(dY,"dY",namedArgs);
	argDouble(dZ,"dZ",namedArgs);
	argDouble(X0,"X0",namedArgs);
	argDouble(Y0,"Y0",namedArgs);
	argDouble(Z0,"Z0",namedArgs);
	argDouble(tolerance,"tolerance",namedArgs);
}

GridImpl::~GridImpl()
{
	if(mapBx) delete mapBx;
	if(mapBy) delete mapBy;
	if(mapBz) delete mapBz;
	if(mapEx) delete mapEx;
	if(mapEy) delete mapEy;
	if(mapEz) delete mapEz;
}

void GridImpl::getFieldValue(const G4double local[4], G4double field[6]) 
									const
{
	G4double x = local[0];
	G4double y = local[1];
	G4double z = local[2];

	x -= X0;
	y -= Y0;
	z -= Z0;

	G4double factor[6];
	factor[0]=factor[1]=factor[2]=factor[3]=factor[4]=factor[5]=1.0;
	if(extendX && x < 0.0) {
		x = -x;
		for(int i=0; i<6; ++i) {
			if(extendXbits & (1<<i)) factor[i] = -factor[i];
		}
	}
	if(extendY && y < 0.0) {
		y = -y;
		for(int i=0; i<6; ++i) {
			if(extendYbits & (1<<i)) factor[i] = -factor[i];
		}
	}
	if(extendZ && z < 0.0) {
		z = -z;
		for(int i=0; i<6; ++i) {
			if(extendZbits & (1<<i)) factor[i] = -factor[i];
		}
	}

	// We compute a 3D linear average of the 8 surrounding points in the map
	// First, get the X,Y,Z indices into i,j,k
	int i = (int)floor(x/dX);
	int j = (int)floor(y/dY);
	int k = (int)floor(z/dZ);
	if(i < 0 || i >= nX-1 || j < 0 || j >= nY-1 || k < 0 || k >= nZ-1) {
		field[0] = field[1] = field[2] = field[3] = field[4] = 
								field[5] = 0.0;
		return;
	}
	// m is the initial index (corner of the cube with minimum X, Y, and Z)
	int m = k*nY*nX + j*nX + i;
	BLAssert(m+nY*nX+nX+1 < nX*nY*nZ);

	// now compute the fractional weighting factors for X, Y, and Z
	float fx = 1.0 - (x - i*dX) / dX;
	BLAssert(fx >= 0.0 && fx <= 1.0);
	float fy = 1.0 - (y - j*dY) / dY;
	BLAssert(fy >= 0.0 && fy <= 1.0);
	float fz = 1.0 - (z - k*dZ) / dZ;
	BLAssert(fz >= 0.0 && fz <= 1.0);

	// now compute the fractional weighting factors for the 8 corners
	float f0 = fx*fy*fz;
	float f1 = (1.0-fx)*fy*fz;
	float f2 = fx*(1.0-fy)*fz;
	float f3 = (1.0-fx)*(1.0-fy)*fz;
	float f4 = fx*fy*(1.0-fz);
	float f5 = (1.0-fx)*fy*(1.0-fz);
	float f6 = fx*(1.0-fy)*(1.0-fz);
	float f7 = (1.0-fx)*(1.0-fy)*(1.0-fz);

	// Finally, compute the components of the field
#define COMPONENT(C)						\
	G4double C = 0.0;					\
	if(map##C) C =	map##C[m]*f0 + map##C[m+1]*f1 + 	\
		map##C[m+nX]*f2 + map##C[m+nX+1]*f3 + 		\
		map##C[m+nY*nX]*f4 + map##C[m+nY*nX+1]*f5 +	\
		map##C[m+nY*nX+nX]*f6 + map##C[m+nY*nX+nX+1]*f7;
	COMPONENT(Bx);
	COMPONENT(By);
	COMPONENT(Bz);
	COMPONENT(Ex);
	COMPONENT(Ey);
	COMPONENT(Ez);

	field[0] = Bx * factor[0];
	field[1] = By * factor[1];
	field[2] = Bz * factor[2];
	field[3] = Ex * factor[3];
	field[4] = Ey * factor[4];
	field[5] = Ez * factor[5];
}

bool GridImpl::handleCommand(InputFile &in, BLArgumentVector &argv, 
			BLArgumentMap &namedArgs)
{
	if(argv[0] == "extendX") {
		extendX = true;
		extendXbits = bits(namedArgs["flip"]);
		return true;
	} else if(argv[0] == "extendY") {
		extendY = true;
		extendYbits = bits(namedArgs["flip"]);
		return true;
	} else if(argv[0] == "extendZ") {
		extendZ = true;
		extendZbits = bits(namedArgs["flip"]);
		return true;
	} else if(argv[0] == "Bx") {
		BLCommand::printError("BLFieldMap: Bx not implemented");
		return true;
	} else if(argv[0] == "By") {
		BLCommand::printError("BLFieldMap: By not implemented");
		return true;
	} else if(argv[0] == "Bz") {
		BLCommand::printError("BLFieldMap: Bz not implemented");
		return true;
	} else if(argv[0] == "Ex") {
		BLCommand::printError("BLFieldMap: Ex not implemented");
		return true;
	} else if(argv[0] == "Ey") {
		BLCommand::printError("BLFieldMap: Ey not implemented");
		return true;
	} else if(argv[0] == "Ez") {
		BLCommand::printError("BLFieldMap: Ez not implemented");
		return true;
	} else if(argv[0] == "data") {
		for(char *line=0; (line=in.getline())!=0; ) {
			if(isalpha(line[0])) {
				in.repeatLine();
				break;
			}
			int n;
			float X=0.0,Y=0.0,Z=0.0,Bx=0.0,By=0.0,Bz=0.0,
							Ex=0.0,Ey=0.0,Ez=0.0;
			for(char *p=line; (p=strchr(p,','))!=0;) *p = ' ';
			n = sscanf(line,"%f%f%f%f%f%f%f%f%f",&X,&Y,&Z,
						&Bx,&By,&Bz,&Ex,&Ey,&Ez);
			if(n <= 3) {
				BLCommand::printError("BLFieldMap: invalid line %d\n",in.linenumber());
				continue;
			}
			setField(X,Y,Z,Bx*tesla,By*tesla,Bz*tesla,
				Ex*megavolt/meter,Ey*megavolt/meter,
				Ez*megavolt/meter,in.linenumber());
		}
		return true;
	} else {
		return false;
	}
}

void GridImpl::getBoundingPoint(int i, G4double point[4])
{
	if(extendX)
		point[0] = (i&1 ? -(nX-1)*dX : (nX-1)*dX) + X0;
	else
		point[0] = (i&1 ? 0.0 : (nX-1)*dX) + X0;
	if(extendY)
		point[1] = (i&2 ? -(nY-1)*dY : (nY-1)*dY) + Y0;
	else
		point[1] = (i&2 ? 0.0 : (nY-1)*dY) + Y0;
	if(extendZ)
		point[2] = (i&4 ? -(nZ-1)*dZ : (nZ-1)*dZ) + Z0;
	else
		point[2] = (i&4 ? 0.0 : (nZ-1)*dZ) + Z0;
}

int GridImpl::bits(G4String s)
{
	int v=0;
	if(s.find("Bx") < s.size()) v |= 1;
	if(s.find("By") < s.size()) v |= 2;
	if(s.find("Bz") < s.size()) v |= 4;
	if(s.find("Ex") < s.size()) v |= 8;
	if(s.find("Ey") < s.size()) v |= 16;
	if(s.find("Ez") < s.size()) v |= 32;
	return v;
}

const char *GridImpl::bits2str(int v)
{
	static char retval[32];
	retval[0] = '\0';
	if(v & 1) strcat(retval,"Bx,");
	if(v & 2) strcat(retval,"By,");
	if(v & 4) strcat(retval,"Bz,");
	if(v & 8) strcat(retval,"Ex,");
	if(v & 16) strcat(retval,"Ey,");
	if(v & 32) strcat(retval,"Ez,");
	return retval;
}

bool GridImpl::setField(G4double X, G4double Y, G4double Z, G4double Bx,
	G4double By, G4double Bz, G4double Ex, G4double Ey, G4double Ez,
	int linenumber)
{
	if(!mapBx && Bx != 0.0) {
		mapBx = new float[nX*nY*nZ];
		BLAssert(mapBx != 0);
		for(int i=0; i<nX*nY*nZ; ++i) 
			mapBx[i] = 0.0;
	}
	if(!mapBy && By != 0.0) {
		mapBy = new float[nX*nY*nZ];
		BLAssert(mapBy != 0);
		for(int i=0; i<nX*nY*nZ; ++i) 
			mapBy[i] = 0.0;
	}
	if(!mapBz && Bz != 0.0) {
		mapBz = new float[nX*nY*nZ];
		BLAssert(mapBz != 0);
		for(int i=0; i<nX*nY*nZ; ++i) 
			mapBz[i] = 0.0;
	}
	if(!mapEx && Ex != 0.0) {
		mapEx = new float[nX*nY*nZ];
		BLAssert(mapEx != 0);
		for(int i=0; i<nX*nY*nZ; ++i) 
			mapEx[i] = 0.0;
	}
	if(!mapEy && Ey != 0.0) {
		mapEy = new float[nX*nY*nZ];
		BLAssert(mapEy != 0);
		for(int i=0; i<nX*nY*nZ; ++i) 
			mapEy[i] = 0.0;
	}
	if(!mapEz && Ez != 0.0) {
		mapEz = new float[nX*nY*nZ];
		BLAssert(mapEz != 0);
		for(int i=0; i<nX*nY*nZ; ++i) 
			mapEz[i] = 0.0;
	}
	X -= X0;
	Y -= Y0;
	Z -= Z0;
	int i = (int)floor((X/dX) + 0.5);
	if(i<0 || fabs(i*dX-X)>tolerance || i >= nX) {
		BLCommand::printError("BLFieldMap: ERROR point off"
			" grid  X=%.2f line=%d\n",X, linenumber);
		return false;
	}
	int j = (int)floor((Y/dY) + 0.5);
	if(j<0 || fabs(j*dY-Y)>tolerance || j >= nY) {
		BLCommand::printError("BLFieldMap: ERROR point off"
			" grid Y=%.2f line=%d\n",Y, linenumber);
		return false;
	}
	int k = (int)floor((Z/dZ) + 0.5);
	if(k<0 || fabs(k*dZ-Z)>tolerance || k >= nZ) {
		BLCommand::printError("BLFieldMap: ERROR point off"
			" grid X=%.2f line=%d\n",Z, linenumber);
		return false;
	}
	int m = k*nY*nX + j*nX + i;
	BLAssert(m >= 0 && m < nX*nY*nZ);
	if(mapBx) mapBx[m] = Bx;
	if(mapBy) mapBy[m] = By;
	if(mapBz) mapBz[m] = Bz;
	if(mapEx) mapEx[m] = Ex;
	if(mapEy) mapEy[m] = Ey;
	if(mapEz) mapEz[m] = Ez;

	return true;
}

bool GridImpl::writeFile(FILE *f)
{
	fprintf(f,"grid X0=%g Y0=%g Z0=%g nX=%d nY=%d nZ=%d dX=%g dY=%g dZ=%g\n",
		X0,Y0,Z0,nX,nY,nZ,dX,dY,dZ);
	if(extendX) fprintf(f,"extendX flip=%s\n",bits2str(extendXbits));
	if(extendY) fprintf(f,"extendY flip=%s\n",bits2str(extendYbits));
	if(extendZ) fprintf(f,"extendZ flip=%s\n",bits2str(extendZbits));
	fprintf(f,"data\n");

	for(int i=0; i<nX; ++i) {
		G4double X = X0 + i*dX;
		for(int j=0; j<nY; ++j) {
			G4double Y = Y0 + j*dY;
			for(int k=0; k<nZ; ++k) {
				G4double Z = Z0 + k*dZ;
				int m = k*nY*nX + j*nX + i;
				BLAssert(m >= 0 && m < nX*nY*nZ);
				G4double Bx = (mapBx ? mapBx[m] : 0.0);
				G4double By = (mapBy ? mapBy[m] : 0.0);
				G4double Bz = (mapBz ? mapBz[m] : 0.0);
				fprintf(f,"%.1f %.1f %.1f %g %g %g",
					X,Y,Z,Bx/tesla,By/tesla,Bz/tesla);
				if(hasE()) {
					G4double Ex = (mapEx ? mapEx[m] : 0.0);
					G4double Ey = (mapEy ? mapEy[m] : 0.0);
					G4double Ez = (mapEz ? mapEz[m] : 0.0);
					fprintf(f," %g %g %g",
						Ex/(megavolt/meter),
						Ey/(megavolt/meter),
						Ez/(megavolt/meter));
				}
				fprintf(f,"\n");
			}
		}
	}
	return true;
}





CylinderImpl::CylinderImpl(BLArgumentVector &argv, BLArgumentMap &namedArgs)
							: FieldMapImpl()
{
	nR = 2;
	nZ = 2;
	dR = 10.0*mm;
	dZ = 10.0*mm;
	Z0 = 0.0;
	tolerance = 0.01*mm;
	mapBz = 0;
	mapBr = 0;
	mapBphi = 0;
	mapEr = 0;
	mapEz = 0;
	extendZ = false;
	extendBrFactor = extendBzFactor = 1.0;
	extendErFactor = extendEzFactor = 1.0;
	argInt(nR,"nR",namedArgs);
	argInt(nZ,"nZ",namedArgs);
	argDouble(dR,"dR",namedArgs);
	argDouble(dZ,"dZ",namedArgs);
	argDouble(Z0,"Z0",namedArgs);
	argDouble(tolerance,"tolerance",namedArgs);
}

CylinderImpl::~CylinderImpl()
{
	if(mapBz) delete mapBz;
	if(mapBr) delete mapBr;
	if(mapBphi) delete mapBphi;
	if(mapEr) delete mapEr;
	if(mapEz) delete mapEz;
}

void CylinderImpl::getFieldValue(const G4double local[4], G4double field[6]) 
									const
{
	G4double z = local[2];
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);

	bool extending = false;
	if(z < 0.0 && extendZ) {
		z = -z;
		extending = true;
	}

	z -= Z0;

	// 2D linear average of the 4 surrounding points in the map
	int i = (int)floor(r/dR);
	int j = (int)floor(z/dZ);
	if(z < 0.0 || i >= nR-1 || j < 0 || j >= nZ-1) {
		field[0] = field[1] = field[2] = field[3] = field[4] = 
								field[5] = 0.0;
		return;
	}

	float fr = 1.0 - (r - i*dR) / dR;
	BLAssert(fr >= 0.0 && fr <= 1.0);
	float fz = 1.0 - (z - j*dZ) / dZ;
	BLAssert(fz >= 0.0 && fz <= 1.0);

	G4double Bz = 0.0;
	if(mapBz) Bz =
		mapBz[j*nR+i]*fr*fz +
		mapBz[j*nR+i+1]*(1.0-fr)*fz +
		mapBz[(j+1)*nR+i]*fr*(1.0-fz) +
		mapBz[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
	G4double Br = 0.0;
	if(mapBr) Br =
		mapBr[j*nR+i]*fr*fz +
		mapBr[j*nR+i+1]*(1.0-fr)*fz +
		mapBr[(j+1)*nR+i]*fr*(1.0-fz) +
		mapBr[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
	G4double Bphi = 0.0;
	if(mapBphi) Bphi =
		mapBphi[j*nR+i]*fr*fz +
		mapBphi[j*nR+i+1]*(1.0-fr)*fz +
		mapBphi[(j+1)*nR+i]*fr*(1.0-fz) +
		mapBphi[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
	G4double Ez = 0.0;
	if(mapEz) Ez =
		mapEz[j*nR+i]*fr*fz +
		mapEz[j*nR+i+1]*(1.0-fr)*fz +
		mapEz[(j+1)*nR+i]*fr*(1.0-fz) +
		mapEz[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
	G4double Er = 0.0;
	if(mapEr) Er =
		mapEr[j*nR+i]*fr*fz +
		mapEr[j*nR+i+1]*(1.0-fr)*fz +
		mapEr[(j+1)*nR+i]*fr*(1.0-fz) +
		mapEr[(j+1)*nR+i+1]*(1.0-fr)*(1.0-fz);
	if(extending) {
		Br *= extendBrFactor;
		Bz *= extendBzFactor;
		Er *= extendErFactor;
		Ez *= extendEzFactor;
	}

	G4double phi = atan2(local[1], local[0]);
	field[0] = Br * cos(phi) - Bphi * sin(phi);
	field[1] = Br * sin(phi) + Bphi * cos(phi);
	field[2] = Bz;
	field[3] = Er * cos(phi);
	field[4] = Er * sin(phi);
	field[5] = Ez;
}

bool CylinderImpl::handleCommand(InputFile &in, BLArgumentVector &argv, 
			BLArgumentMap &namedArgs)
{
	if(argv[0] == "extendZ") {
		extendZ = true;
		const char *p = namedArgs["flip"].c_str();
		if(strstr(p,"Br")) extendBrFactor = -1.0;
		if(strstr(p,"Bz")) extendBzFactor = -1.0;
		if(strstr(p,"Er")) extendErFactor = -1.0;
		if(strstr(p,"Ez")) extendEzFactor = -1.0;
		return true;
	} else if(argv[0] == "Br") {
		if(mapBr) return false;
		mapBr = new float[nR*nZ];
		for(int i=0; i<nR*nZ; ++i) mapBr[i] = 0.0;
		return readBlock(in,mapBr,nZ,nR,tesla);
	} else if(argv[0] == "Bz") {
		if(mapBz) return false;
		mapBz = new float[nR*nZ];
		for(int i=0; i<nR*nZ; ++i) mapBz[i] = 0.0;
		return readBlock(in,mapBz,nZ,nR,tesla);
	} else if(argv[0] == "Er") {
		if(mapEr) return false;
		mapEr = new float[nR*nZ];
		for(int i=0; i<nR*nZ; ++i) mapEr[i] = 0.0;
		return readBlock(in,mapEr,nZ,nR,megavolt/meter);
	} else if(argv[0] == "Ez") {
		if(mapEz) return false;
		mapEz = new float[nR*nZ];
		for(int i=0; i<nR*nZ; ++i) mapEz[i] = 0.0;
		return readBlock(in,mapEz,nZ,nR,megavolt/meter);
	} else if(argv[0] == "data") {
		for(char *line=0; (line=in.getline())!=0;) {
			if(isalpha(line[0])) {
				in.repeatLine();
				break;
			}
			int n;
			float R=0.0,Z=0.0,Br=0.0,Bz=0.0,Er=0.0,Ez=0.0;
			for(char *p=line; (p=strchr(p,','))!=0;) *p = ' ';
			n = sscanf(line,"%f%f%f%f%f%f",&R,&Z,&Br,&Bz,&Er,&Ez);
			if(n <= 2) continue;
			setField(R,Z,Bz*tesla,Br*tesla,0.0,Er*megavolt/meter,
				Ez*megavolt/meter,in.linenumber());
		}
		return true;
	} else {
		return false;
	}
}

void CylinderImpl::getBoundingPoint(int i, G4double point[4])
{
	point[0] = (i&1 ? 1.0 : -1.0) * (nR-1) * dR;
	point[1] = (i&2 ? 1.0 : -1.0) * (nR-1) * dR;
	if(extendZ)
		point[2] = (i&4 ? 1.0 : -1.0) * (nZ-1) * dZ;
	else
		point[2] = (i&4 ? Z0 : Z0+(nZ-1)*dZ);
}

bool CylinderImpl::setField(G4double R, G4double Z, G4double Bz, G4double Br,
		G4double Bphi, G4double Er, G4double Ez, int linenumber)
{
	if((Br != 0.0 || Bz != 0.0 || Bphi != 0.0) && !mapBr) {
		mapBr = new float[nR*nZ];
		mapBz = new float[nR*nZ];
		mapBphi = new float[nR*nZ];
		BLAssert(mapBr != 0 && mapBz != 0 && mapBphi != 0);
		for(int i=0; i<nR*nZ; ++i) 
			mapBz[i] = mapBr[i] = mapBphi[i] = 0.0;
	}
	if((Er != 0.0 || Ez != 0.0) && !mapEr) {
		mapEr = new float[nR*nZ];
		mapEz = new float[nR*nZ];
		BLAssert(mapEr != 0 && mapEz != 0);
		for(int i=0; i<nR*nZ; ++i) 
			mapEr[i] = mapEz[i] = 0.0;
	}
	int i = (int)floor((R/dR) + 0.5);
	if(i<0 || fabs(i*dR-R)>tolerance || i >= nR) {
		BLCommand::printError("BLFieldMap: ERROR point off"
			    " grid R=%.1f line=%d\n",R,linenumber);
		return false;
	}
	int j = (int)floor(((Z-Z0)/dZ) + 0.5);
	if(j<0 || fabs(j*dZ+Z0-Z)>tolerance || j >= nZ) {
		BLCommand::printError("BLFieldMap: ERROR point off"
			    " grid Z=%.1f line=%d\n",Z,linenumber);
		return false;
	}
	if(mapBr) {
		mapBz[j*nR+i] = Bz;
		mapBr[j*nR+i] = Br;
		mapBphi[j*nR+i] = Bphi;
	}
	if(mapEr) {
		mapEr[j*nR+i] = Er;
		mapEz[j*nR+i] = Ez;
	}

	return true;
}

bool CylinderImpl::writeFile(FILE *f)
{
	fprintf(f,"cylinder Z0=%g nR=%d nZ=%d dR=%g dZ=%g\n",Z0,nR,nZ,dR,dZ);
	if(extendZ) {
		fprintf(f,"extendZ flip=");
		if(extendBrFactor < 0.0) fprintf(f,"Br,");
		if(extendBzFactor < 0.0) fprintf(f,"Bz,");
		if(extendErFactor < 0.0) fprintf(f,"Er,");
		if(extendEzFactor < 0.0) fprintf(f,"Ez,");
		fprintf(f,"\n");
	}
	fprintf(f,"data\n");

	for(int i=0; i<nR; ++i) {
		G4double R = i*dR;
		for(int j=0; j<nZ; ++j) {
			G4double Z = Z0 + j*dZ;
			int m = j*nR + i;
			BLAssert(m >= 0 && m < nR*nZ);
			G4double Br = (mapBr ? mapBr[m] : 0.0);
			G4double Bz = (mapBz ? mapBz[m] : 0.0);
			fprintf(f,"%.1f %.1f %g %g", R,Z,Br/tesla,Bz/tesla);
			if(hasE()) {
				G4double Er = (mapEr ? mapEr[m] : 0.0);
				G4double Ez = (mapEz ? mapEz[m] : 0.0);
				fprintf(f," %g %g", Er/(megavolt/meter),
						    Ez/(megavolt/meter));
			}
			fprintf(f,"\n");
		}
	}
	return true;
}


