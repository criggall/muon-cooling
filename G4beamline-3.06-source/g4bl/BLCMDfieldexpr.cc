//	BLCMDfieldexpr.cc
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
#include "BLCommand.hh"
#include "BLEvaluator.hh"
#include "BLTune.hh"
#include "BLElement.hh"
#include "BLGlobalField.hh"
#include "BLElementField.hh"

/**	class BLCMDfieldexpr implements a fieldexpr command
 *
 *	This command has no volumes, logical or physical. it just implements
 *	a field (E and/or B) from expressions on its command line.
 **/
class BLCMDfieldexpr : public BLElement {
	G4double factorB;
	G4double factorE;
	G4double timeOffset;
	G4String Bx;
	G4String By;
	G4String Bz;
	G4String Br;
	G4String Bphi;
	G4String Ex;
	G4String Ey;
	G4String Ez;
	G4String Er;
	G4String time;
	G4double length;
	G4double width;
	G4double height;
	G4double radius;
	G4double tmin;
	G4double tmax;
	int nX;
	int nY;
	int nZ;
	int nR;
	int nT;
	G4double tolerance;
	BLFieldMap *map;
	friend class FieldComputation;
	void handleTimeDependence();
	BLEvaluator eval;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDfieldexpr();

	/// Destructor.
	virtual ~BLCMDfieldexpr() { }

	/// Copy constructor.
	BLCMDfieldexpr(const BLCMDfieldexpr& r);

	/// clone()
	BLElement *clone() { return new BLCMDfieldexpr(*this); }

	/// commandName() returns "fieldexpr".
	G4String commandName() { return "fieldexpr"; }

	/// command() implements the fieldexpr command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the fieldexpr
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns length (so place front=1 will work)
	G4double getLength() { return length; }

	/// getWidth() returns 0 (no physical volume nere)
	G4double getWidth() { return 0.0; }

	/// getHeight() returns 0 (no physical volume nere)
	G4double getHeight() { return 0.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement. (no volume => every point is "outside")
	bool isOutside(G4ThreeVector &local, G4double tolerance) 
		{ return true; }

	/// generatePoints() from BLElement. (no volume => no generate)
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v) 
		{ v.clear(); }

	/// maxError() returns the maximum relative error 
	/// (i.e. error / max field).
	G4double maxError(class FieldComputation *p);

	/// evaluateAndCheck() evaluates an expression and checks it for 
	/// validity. Returns 0.0 for an empty string.
	G4double evaluateAndCheck(G4String e) {
		if(e.size() == 0) return 0.0;
		G4double v=eval.evaluate(e.c_str());
		if(!eval.isOK()) {
			G4Exception("fieldexpr","Invalid Expression",
				FatalException,e.c_str());
		}
		return v;
	}

	/// setVariable() sets a variable value for evaluateAndCheck()
	void setVariable(const char *_name, G4double v) {
		eval.setVariable(_name,v);
	}

	/// checkValidExpr() checks that the evaluator had a valid expression.
	/// Prints error message and exits the program if not.
	void checkValidExpr(BLEvaluator &e, G4String s)
		{ if(!e.isOK()) {
			G4Exception("fieldexpr","Invalid Expression",
				FatalException,s.c_str());
		  }
		}
};

BLCMDfieldexpr defaultBLCMDfieldexpr;	// default object

class FieldExprPlacement : public BLElementField {
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double *factorB;
	G4double *factorE;
	G4double timeOffset;
	BLFieldMap *map;
public:
	FieldExprPlacement(BLFieldMap *_map, G4RotationMatrix *rot,
		G4ThreeVector &pos, G4double& _factorB, G4double& _factorE,
		G4double _timeOffset);
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};

class FieldComputation : public G4ElectroMagneticField {
	BLCMDfieldexpr *expr;
public:
	FieldComputation(BLCMDfieldexpr *_expr) : G4ElectroMagneticField() 
		{ expr=_expr; }
	G4bool DoesFieldChangeEnergy() const { return true; }
	void GetFieldValue(const G4double point[4], G4double *field) const;
};

// Default constructor - be sure to use the default constructor BLElement()
BLCMDfieldexpr::BLCMDfieldexpr() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("implements a field map, E and/or B, from expressions.");
	setDescription("A fieldexpr element can be either a box or a cylinder; "
		"set length and radius for cylinder, set length and width and "
		"height for a box. Units are Tesla, MegaVolts/meter, mm, and "
		"ns. Expressions for the field components can use {x,y,z} for "
		"a box or {z,r} for a cylinder; the time expression can use "
		"{t}. If present, the time expression multiples all "
		"components.\n\n"
		"Expressions can use all C operators except ?:, and x^n is x "
		"to the nth power (n integer). The following functions are "
		"available: abs(), min(),max(), sqrt(), pow(), sin(), cos(), "
		"tan(), asin(), acos(), atan(), atan2(), sinh(), cosh(), "
		"tanh(), exp(), log(), log10(), floor(), ceil(), if(). "
		"if(condition,a,b) replaces the C (condition? a : b).\n\n"
		"A field map is used for tracking efficiency; the "
		"number of points in the map is increased "
		"until the largest map error divided by the maximum field "
		"is smaller than tolerance, or 1M points is exceeded. "
		"Similarly for the time dependence.\n\n"
		"For time dependence: if t-timeOffset<tmin, the value at tmin "
		"is used; if t-timeOffset>tmax, the value at tmax is used.\n\n"
		"Note that divide by zero is reported as invalid expression. "
		"For a Li lens you probably want to use an expression like "
		"this: 'if(r<100,500.0*r/100,500.0*100/r)', where 100mm is "
		"the radius, and 500T is the field at r=100mm.\n\n"
		"This field must be placed (via the place command); that "
		"specifies where (x=0,y=0,z=0) of the map is located in the "
		"parent. "
		"Note that front=1 can be used to place this field.");

	// provide initial values for fields
	factorB = 1.0;
	factorE = 1.0;
	timeOffset = 0.0;
	Bx = By = Bz = Br = Bphi = Ex = Ey = Ez = Er = time = "";
	length = width = height = tmin = tmax = 0.0;
	radius = -1.0;
	nX = nY = nZ = nR = nT = 11;
	tolerance = 0.001;
	map = 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDfieldexpr::BLCMDfieldexpr(const BLCMDfieldexpr& r) : BLElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	BLTune::copyTunableArg(&factorB,&r.factorB);
	BLTune::copyTunableArg(&factorE,&r.factorE);
	timeOffset = r.timeOffset;
	Bx = r.Bx;
	By = r.By;
	Bz = r.Bz;
	Br = r.Br;
	Bphi = r.Bphi;
	Ex = r.Ex;
	Ey = r.Ey;
	Ez = r.Ez;
	Er = r.Er;
	time = r.time;
	length = r.length;
	width = r.width;
	height = r.height;
	radius = r.radius;
	tmin = r.tmin;
	tmax = r.tmax;
	nX = r.nX;
	nY = r.nY;
	nZ = r.nZ;
	nR = r.nR;
	nT = r.nT;
	tolerance = r.tolerance;
	map = r.map;
}

int BLCMDfieldexpr::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		BLCommand::printError("fieldexpr: Invalid command, must have one name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultBLCMDfieldexpr.handleNamedArgs(namedArgs);
	}

	BLCMDfieldexpr *t = new BLCMDfieldexpr(defaultBLCMDfieldexpr);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	FieldComputation *p = new FieldComputation(t);
	t->map = new BLFieldMap();

	// check error; if necessary, double # map points in each dimension.
	double err=0.0;
	int ntot=0;
	for(;;) {
	    if(t->radius > 0.0) {
		t->map->createCylinderMap(-t->length/2,t->radius/(t->nR-1),
			t->length/(t->nZ-1),t->nR,t->nZ,p);
		ntot = t->nR*t->nZ;
	    } else {
		t->map->createGridMap(-t->width/2,-t->height/2,-t->length/2,
			t->width/(t->nX-1),t->height/(t->nY-1),
			t->length/(t->nZ-1),t->nX,t->nY,t->nZ,p);
		ntot = t->nX*t->nY*t->nZ;
	    }
	    err = t->maxError(p);
	    printf("fieldexpr %s: Map %d points  Max Relative Error = %.4f\n",
			    t->getName().c_str(),ntot,err);
	    if(err <= t->tolerance || ntot >= 1000000)
		    break;
	    if(t->radius > 0.0) {
		t->nR = (t->nR-1)*2 + 1;
		t->nZ = (t->nZ-1)*2 + 1;
	    } else {
		t->nX = (t->nX-1)*2 + 1;
		t->nY = (t->nY-1)*2 + 1;
		t->nZ = (t->nZ-1)*2 + 1;
	    }
	}

	delete p;

	t->handleTimeDependence();

	t->print(argv[0]);

	return retval;
}

void BLCMDfieldexpr::defineNamedArgs()
{
	argTunable(factorB,"factorB","Factor for the B-field (1.0).");
	argTunable(factorE,"factorE","Factor for the E-field (1.0).");
	argDouble(timeOffset,"timeOffset","Time offset (ns).");
	// the rest are not permitted to change
	argString(Bx,"Bx","Expression for Bx (Tesla), use {x,y,z}.",false);
	argString(By,"By","Expression for By (Tesla), use {x,y,z}.",false);
	argString(Bz,"Bz","Expression for Bz (Tesla), use {x,y,z} or {r,z}.",false);
	argString(Br,"Br","Expression for Br (Tesla), use {r,z}.",false);
	argString(Bphi,"Bphi","Expression for Bphi (Tesla), use {r,z}.",false);
	argString(Ex,"Ex","Expression for Ex (MV/m), use {x,y,z}.",false);
	argString(Ey,"Ey","Expression for Ey (MV/m), use {x,y,z}.",false);
	argString(Ez,"Ez","Expression for Ez (MV/m), use {x,y,z} or {r,z}.",false);
	argString(Er,"Er","Expression for Er (MV/m), use {r,z}.",false);
	argString(time,"time","Expression for time-dependence factor, use {t}.",false);
	argInt(nX,"nX","Number of grid points in x.",false);
	argInt(nY,"nY","Number of grid points in y.",false);
	argInt(nZ,"nZ","Number of grid points in z.",false);
	argInt(nR,"nR","Number of grid points in r.",false);
	argInt(nT,"nT","Number of grid points in t.",false);
	argDouble(tolerance,"tolerance","Required relative accuracy (0.001).",1,"",false);
	argDouble(length,"length","Length of field map (mm).",mm,"",false);
	argDouble(width,"width","Width of rectangular field map (mm).",mm,"",false);
	argDouble(height,"height","Height of rectangular field map (mm).",mm,"",false);
	argDouble(radius,"radius","Radius of cylindrical field map (mm).",mm,"",false);
	argDouble(tmin,"tmin","Minimum value of t (ns).",mm,"",false);
	argDouble(tmax,"tmax","Maximum value of t (ns).",mm,"",false);
}

void BLCMDfieldexpr::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	// get globalRotation and globalPosition
	G4RotationMatrix *globalRotation = 0;
	if(relativeRotation && parentRotation) {
		globalRotation = 
		    new G4RotationMatrix(*parentRotation * *relativeRotation);
	} else if(relativeRotation) {
		globalRotation = relativeRotation;
	} else if(parentRotation) {
		globalRotation = parentRotation;
	}
	G4ThreeVector globalPosition(relativePosition + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	FieldExprPlacement *p = new FieldExprPlacement(map,globalRotation,
				globalPosition,factorB,factorE,timeOffset);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLMappedMagnet::Construct %s parent=%s relZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2]);
}

void BLCMDfieldexpr::handleTimeDependence()
{
	if(time != "" && nT > 1) {
		static BLEvaluator e;
		double *t, *f;
		for(;;) {
			double dt = (tmax - tmin) / (nT-1);
			t = new double[nT];
			f = new double[nT];
			for(int i=0; i<nT; ++i) {
				t[i] = tmin + dt*i;
				e.setVariable("t",t[i]);
				f[i] = e.evaluate(time);
				checkValidExpr(e,time+"A");
			}
			// check error, double nT if needed
			double err=0.0;
			double v=tmin+dt/2.0;
			for(int i=0; i<nT-1; ++i) {
				e.setVariable("t",v);
				double d = (f[i]+f[i+1])/2.0 - e.evaluate(time);
				checkValidExpr(e,time+"B");
				if(fabs(d) > err) err = fabs(d);
				v += dt;
			}
	    		printf("fieldexpr %s: Time %d points  Max Relative Error=%.4f\n",
			    getName().c_str(),nT,err);
			if(err <= tolerance || nT > 1000000)
				break;
			nT = (nT-1)*2 + 1;
			delete[] t;
			delete[] f;
		}
		map->createTimeDependence(nT,t,f,f);
		delete[] t;
		delete[] f;
	}
}

G4double BLCMDfieldexpr::maxError(FieldComputation *p)
{
	G4double err[6],maxField[6];
	for(int i=0; i<6; ++i)
		err[i] = maxField[i] = 0.0;

	// find the max error at the center of every grid box
	if(radius > 0.0) {
		double dR = radius/(nR-1);
		double dZ = length/(nZ-1);
		for(double z=(-length+dZ)/2.0; z<length/2.0; z+=dZ) {
			for(double r=dR/2.0; r<radius; r+=dR) {
				G4double point[4],field1[6],field2[6];
				point[0] = r;
				point[1] = 0.0;
				point[2] = z;
				point[3] = 0.0; // map has no time dependence
				map->getFieldValue(point,field1,1.0,1.0);
				p->GetFieldValue(point,field2);
				for(int i=0; i<6; ++i) {
					if(fabs(field2[i]) > maxField[i])
					    maxField[i] = fabs(field2[i]);
					if(fabs(field2[i]-field1[i]) > err[i])
					    err[i] = fabs(field2[i]-field1[i]);
				}
			}
		}
	} else {
		double dX = width/(nX-1);
		double dY = height/(nY-1);
		double dZ = length/(nZ-1);
		for(double z=(-length+dZ)/2.0; z<length/2.0; z+=dZ) {
		    for(double y=(-height+dY)/2.0; y<height/2.0; y+= dY) {
			for(double x=(-width+dX)/2.0; x<width/2.0; x+=dX) {
			    G4double point[4],field1[6],field2[6];
			    point[0] = x;
			    point[1] = y;
			    point[2] = z;
			    point[3] = 0.0; // map has no time dependence, yet
			    map->getFieldValue(point,field1,1.0,1.0);
			    p->GetFieldValue(point,field2);
//@ printf("fieldexpr x,y,z=%.3f,%.3f,%.3f   B1=%.4f,%.4f,%.4f   B2=%.4f,%.4f,%.4f\n",x,y,z,field1[0]/tesla,field1[1]/tesla,field1[2]/tesla,field2[0]/tesla,field2[1]/tesla,field2[2]/tesla);

			    for(int i=0; i<6; ++i) {
				if(fabs(field2[i]) > maxField[i])
					maxField[i] = fabs(field2[i]);
				if(fabs(field2[i]-field1[i]) > err[i])
					err[i] = fabs(field2[i]-field1[i]);
			    }
			}
		    }
		}
	}

	// denominator is largest component, not largest norm -- good enough.
	if(maxField[1] > maxField[0]) maxField[0] = maxField[1];
	if(maxField[2] > maxField[0]) maxField[0] = maxField[2];
	if(maxField[4] > maxField[3]) maxField[3] = maxField[4];
	if(maxField[5] > maxField[3]) maxField[3] = maxField[5];
	if(err[1] > err[0]) err[0] = err[1];
	if(err[2] > err[0]) err[0] = err[2];
	if(err[4] > err[3]) err[3] = err[4];
	if(err[5] > err[3]) err[3] = err[5];
	G4double errB=0.0, errE=0.0;
	if(maxField[0] > 0.000001) errB = err[0]/maxField[0];
	if(maxField[3] > 0.000001) errE = err[3]/maxField[3];
	return (errB>errE ? errB : errE);
}


FieldExprPlacement::FieldExprPlacement(BLFieldMap *_map, G4RotationMatrix *rot,
		G4ThreeVector &pos, G4double& _factorB, G4double& _factorE,
		G4double _timeOffset)
			: global2local(rot,pos), rotation() 
{
	factorB = &_factorB;
	factorE = &_factorE;
	timeOffset = _timeOffset;
	map = _map;

	if(global2local.isRotated()) {
		rotation = global2local.getRotation();
		rotation.invert();
	}

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	for(int i=0; i<8; ++i) {
		map->getBoundingPoint(i,local);
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
	}
}

void FieldExprPlacement::addFieldValue(const G4double point[4],
						G4double field[6]) const
{
	G4double local[4], thisField[6];

	global2local.getLocal(local,point);

	local[3] -= timeOffset;

	map->getFieldValue(local,thisField,*factorB,*factorE);

	if(map->hasB()) {
		if(!rotation.isIdentity()) {
			G4ThreeVector B(thisField[0],thisField[1],thisField[2]);
			B = rotation * B;
			field[0] += B[0];
			field[1] += B[1];
			field[2] += B[2];
		} else {
			field[0] += thisField[0];
			field[1] += thisField[1];
			field[2] += thisField[2];
		}
	}

	if(map->hasE()) {
		if(!rotation.isIdentity()) {
			G4ThreeVector E(thisField[3],thisField[4],thisField[5]);
			E = rotation * E;
			field[3] += E[0];
			field[4] += E[1];
			field[5] += E[2];
		} else {
			field[3] += thisField[3];
			field[4] += thisField[4];
			field[5] += thisField[5];
		}
	}
}

void FieldComputation::GetFieldValue(const G4double point[4], G4double *field) const
{
	if(expr->radius > 0.0) {
		// using X-Z plane (i.e. Y=0)
		BLAssert(fabs(point[1]) < 0.000001 && point[0] >= 0.0);
		expr->setVariable("r",point[0]);
		expr->setVariable("z",point[2]);
		field[0] = expr->evaluateAndCheck(expr->Br) * tesla;
		field[1] = expr->evaluateAndCheck(expr->Bphi) * tesla;
		field[3] = expr->evaluateAndCheck(expr->Er) * (megavolt/meter);
		field[4] = 0.0;
	} else {
		expr->setVariable("x",point[0]);
		expr->setVariable("y",point[1]);
		expr->setVariable("z",point[2]);
		field[0] = expr->evaluateAndCheck(expr->Bx) * tesla;
		field[1] = expr->evaluateAndCheck(expr->By) * tesla;
		field[3] = expr->evaluateAndCheck(expr->Ex) * (megavolt/meter);
		field[4] = expr->evaluateAndCheck(expr->Ey) * (megavolt/meter);
	}
	field[2] = expr->evaluateAndCheck(expr->Bz) * tesla;
	field[5] = expr->evaluateAndCheck(expr->Ez) * (megavolt/meter);

	return;
}
