//	BLEngeFunction.hh
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

#ifndef BLENGEFUNCTION_HH
#define BLENGEFUNCTION_HH

#include <math.h>

/**	enum BLEngeType specifies the default parameters of BLEngeFunction.
 **/
enum BLEngeType { ENGE_BLOCK, ENGE_BEND, ENGE_QUAD, ENGE_OTHER };

/**	class EngeFunction implements the Enge function for fringe fields.
 *
 *	z is the distance from the nominal edge, with z=0 being the edge.
 *	it should be divided by the aperture diameter or full width/height.
 *	z<0 is inside, z>0 is outside.
 *
 *	See the COSY reference manual (pp 32-35) for suggested values of
 *	a1-a6, or use the ENGE_BEND or ENGE_QUAD types (which come from there).
 *	http://cosy.pa.msu.edu/cosymanu/index.html
 *
 *	Mathematica was used to compute the derivatives.
 **/
class BLEngeFunction {
	BLEngeType type;
	double a1,a2,a3,a4,a5,a6;
public:
	/// default constructor.
	BLEngeFunction() { type=ENGE_BLOCK; set(0,0,0,0,0,0); }
	/// constructor for common magnet types.
	BLEngeFunction(BLEngeType t) {
		switch(t) {
		case ENGE_BLOCK:
		case ENGE_OTHER:
		  set(0,0,0,0,0,0);
		  break;
		case ENGE_BEND:
		  set(0.478959,1.911289,-1.185953,1.630554,-1.082657,0.318111);
		  break;
		case ENGE_QUAD:
		  set(0.296471,4.533219,-2.270982,1.068627,-0.036391,0.022261);
		  break;
		}
		type = t;
	}
	/// general constructor.
	BLEngeFunction(double _a1, double _a2, double _a3, double _a4,
						double _a5, double _a6) {
		set(_a1,_a2,_a3,_a4,_a5,_a6);
	}
	/// set the parameters.
	void set(double _a1, double _a2, double _a3, double _a4,
						double _a5, double _a6) {
		a1=_a1; a2=_a2; a3=_a3; a4=_a4; a5=_a5; a6=_a6;
		if(a1==0.0 && a2==0.0 && a3==0.0 && a4==0.0 && a5==0.0 && 
								a6==0.0)
			type = ENGE_BLOCK;
		else
			type = ENGE_OTHER;
	}
	/// evaluate the Enge function at z.
	double operator()(double z) const {
		if(type == ENGE_BLOCK) return (z<=0.0 ? 1.0 : 0.0);
		if(z < -4.0) return 1.0;
		if(z > 4.0) return 0.0;
		return 1.0/(1.0+exp(a1+z*(a2+z*(a3+z*(a4+z*(a5+z*a6))))));
	}
	/// evaluate the derivative of the Enge function at z.
	double prime(double z) const {
		if(type == ENGE_BLOCK) return 0.0;
		if(fabs(z) > 4.0) return 0.0;
		double exp1 = exp(a1+z*(a2+z*(a3+z*(a4+z*(a5+z*a6)))));
		return -exp1/(1.0+exp1)/(1.0+exp1)*
			(a2+z*(2.0*a3+z*(3.0*a4+z*(4.0*a5+z*5.0*a6))));
	}
	double first(double z) { return prime(z); }
	/// evaluate the second derivative of the Enge function at z.
	double second(double z) const {
		if(type == ENGE_BLOCK) return 0.0;
		if(fabs(z) > 4.0) return 0.0;
		double f1 = a1+z*(a2+z*(a3+z*(a4+z*(a5+z*a6))));
		double f2 = (a2+2*a3*z+3*a4*z*z+4*a5*z*z*z+5*a6*z*z*z*z);
		double f3 = (2*a3+6*a4*z+12*a5*z*z+20*a6*z*z*z);
		double exp1 = exp(f1);
		return exp1*((exp1-1.0)*f2*f2-(1.0+exp1)*f3)/
					(1.0+exp1)/(1.0+exp1)/(1.0+exp1);
	}
	/// evaluate the third derivative of the Enge function at z.
	double third(double z) const {
		if(type == ENGE_BLOCK) return 0.0;
		if(fabs(z) > 4.0) return 0.0;
		double f1 = a1+z*(a2+z*(a3+z*(a4+z*(a5+z*a6))));
		double f2 = a2+z*(2*a3+z*(3*a4+4*a5*z+5*a6*z*z));
		double f3 = 2*(a3+z*(3*a4+2*z*(3*a5+5*a6*z)));
		double f4 = a4+2.0*z*(2.0*a5+5.0*a6*z);
		double exp1 = exp(f1);
		double onepexp1 = 1.0 + exp1;
		return -exp1*(6*exp1*exp1*f2*f2*f2-6*exp1*f2*(f2*f2+f3)*onepexp1
			+(f2*f2*f2+3*f2*f3+6*f4)*onepexp1*onepexp1)
			/(onepexp1*onepexp1*onepexp1*onepexp1);
	}
	/// evaluate the fourth derivative of the Enge function at z.
	double fourth(double z) const {
		if(type == ENGE_BLOCK) return 0.0;
		if(fabs(z) > 4.0) return 0.0;
		double f1 = a1+z*(a2+z*(a3+z*(a4+z*(a5+z*a6))));
		double f2 = a2+z*(2*a3+z*(3*a4+4*a5*z+5*a6*z*z));
		double f3 = 2*(a3+z*(3*a4+2*z*(3*a5+5*a6*z)));
		double f4 = a4+2.0*z*(2.0*a5+5.0*a6*z);
		double f5 = a5 + 5*a6*z;
		double exp1 = exp(f1);
		double onepexp1 = 1.0 + exp1;
		return -exp1*(-24*exp1*exp1*exp1*f2*f2*f2*f2+onepexp1*
		    (36*exp1*exp1*f2*f2*(f2*f2+f3)-2*exp1*(7*f2*f2*f2*f2
		    +18*f2*f2*f3+3*f3*f3+24*f2*f4)*onepexp1
		    +(f2*f2*f2*f2+6*f2*f2*f3+3*f3*f3+24*f2*f4+24*f5)
		    *onepexp1*onepexp1))
		    /(onepexp1*onepexp1*onepexp1*onepexp1*onepexp1);
	}
	/// evaluate the fifth derivative of the Enge function at z.
	double fifth(double z) const {
		if(type == ENGE_BLOCK) return 0.0;
		if(fabs(z) > 4.0) return 0.0;
		double f1 = a1+z*(a2+z*(a3+z*(a4+z*(a5+z*a6))));
		double f2 = a2+z*(2*a3+z*(3*a4+4*a5*z+5*a6*z*z));
		double f3 = 2*(a3+z*(3*a4+2*z*(3*a5+5*a6*z)));
		double f4 = a4+2.0*z*(2.0*a5+5.0*a6*z);
		double f5 = a5 + 5*a6*z;
		double exp1 = exp(f1);
		double onepexp1 = 1.0 + exp1;
		return -exp1/(onepexp1*onepexp1*onepexp1*onepexp1*onepexp1*onepexp1)
		    *(120*exp1*exp1*exp1*exp1*f2*f2*f2*f2*f2
		    -240*exp1*exp1*exp1*f2*f2*f2*(f2*f2+f3)*onepexp1
		    +onepexp1*onepexp1*(30*exp1*exp1*f2*(5*f2*f2*f2*f2
		    +12*f2*f2*f3+3*f3*f3+12*f2*f4)-10*exp1*(3*f2*f2*f2*f2*f2
		    +14*f2*f2*f2*f3+9*f2*f3*f3+36*f2*f2*f4+12*f3*f4+24*f2*f5)
		    *onepexp1+(120*a6+f2*f2*f2*f2*f2+10*f2*f2*f2*f3+60*f2*f2*f4
		    +60*f3*f4+15*f2*(f3*f3+8*f5))*onepexp1*onepexp1));
	}
	/// return the type of Enge function
	BLEngeType getType() const { return type; }
};

#endif // BLENGEFUNCTION_HH
