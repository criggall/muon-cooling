//	BLEvaluator.hh
//#define UNARY_MINUS_OK // As of CLHEP 2.1.3.1 it is OK, 2.1.1.0 not

#ifndef BLEVALUATOR_HH
#define BLEVALUATOR_HH

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits>
#include "CLHEP/Evaluator/Evaluator.h"

#ifndef NO_GEANT4	// for testing basic functionality
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "BLCoordinates.hh"
#include "BLGlobalField.hh"
#include "BLAssert.hh"
#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;
#else
#include <assert.h>
#define BLAssert(expr) assert(expr)
#endif

/** class BLEvaluator implements an expression evaluator.
 *
 *	Uses HepTool::Evaluator, but fixes up a few minor bugs:
 *		automatically calls setStdMath()
 *		adds floor(x), ceil(x), if(test,trueVal,falseVal)
 *		fixes unary +/- by enclosing them in parens
 *
 *	NOTE: This is expensive, and should not be used during tracking unless
 *	unavoidable.
 *
 *	This class is entirely inline, no .cc file is used.
 **/
class BLEvaluator : public HepTool::Evaluator {
static double iffunc(double condition, double truevalue, double falsevalue) {
	return condition ? truevalue : falsevalue; }
public:
	/// Constructor
	BLEvaluator() : HepTool::Evaluator() {
		setSystemOfUnits(1.e+3, 1./1.60217733e-25, 1.e+9,
				1./1.60217733e-10, 1.0, 1.0, 1.0);
		setStdMath();
		setFunction("floor",floor);
		setFunction("ceil",ceil);
		setFunction("if",iffunc);
	}

	/// evaluate() will evaluate the expression.
	/// returns NAN if the expression is not valid (use either isOK() or
	/// std::isnan(double) from <cmath> to check).
#ifdef UNARY_MINUS_OK
	double evaluate(const char *e) {
		double v = HepTool::Evaluator::evaluate(e);
		if(status() != OK) v = std::numeric_limits<double>::quiet_NaN();
		return v;
	}
#else //UNARY_MINUS_OK
	double evaluate(const char *expression) {
		char out[1024];
		BLAssert(strlen(expression) < sizeof(out));
		// remove whitespace (simpler parsing below)
		int j=0;
		for(int i=0; expression[i]!='\0'; ++i)
			if(!isspace(expression[i])) out[j++] = expression[i];
		out[j] = '\0';
		// Workaround for serious bug in HepTool::Evaluator --
		// it can only handle a unary +/- after a '(', so
		// put parens around unary +/-
		bool again=false;
		do {
		    again = false;
		    //@ char *in = strdup(out); -- not in Windows
		    char *in = new char[strlen(out)+1];
		    strcpy(in,out); // copies out to in
		    //@ 
		    j = 0;
		    char prev=' ';
		    for(int i=0; in[i]!='\0'; ++i) {
			char c = in[i];
			if((c == '-' || c == '+') && isopchar(prev)) {
				again = true;
				out[j++] = '(';
				out[j++] = c;
				++i;
				int level = 0;
				while(isumchar(in[i]) || level > 0) {
				    if(in[i] == '(')
					++level;
				    if(in[i] == ')')
					--level;
				    // unary +/- after operator -- copy now,
				    // handle during next loop (again = true)
				    // (^ is the only op with higher precedence)
				    if(in[i]=='^' && (in[i+1]=='-' ||
				    				in[i+1]=='+'))
					out[j++] = in[i++];
				    out[j++] = in[i++];
				}
				c = ')';
				--i;
			}
			prev = out[j++] = c;
		    }
		    BLAssert((unsigned)j < sizeof(out));
		    out[j] = '\0';
		    // if empty or all whitespace, force invalid
		    if(j == 0) strcpy(out,"x+-*/y");
		    free(in);
		} while(again);

		double v = HepTool::Evaluator::evaluate(out);
		if(status() != OK) v = std::numeric_limits<double>::quiet_NaN();
//printf("BLEvaluator::evaluate: expr='%s' out='%s' v=%.6f\n",expression,out,v);
		return v;
	}
#endif //UNARY_MINUS_OK

	/// isOK() returns true if the previous call to evaluate() succeeded.
	bool isOK() { return status() == OK; }

#ifndef NO_GEANT4
	/// setTrackVariables() devines variables in the evaluator
	/// corresponding to the usual track variables. These can then
	/// be used in expressions.
	void setTrackVariables(const G4Track *track, BLCoordinateType coordType,
					G4String suffix="", bool fields=false) {
		G4RunManager* runmgr = G4RunManager::GetRunManager();
		const G4Event* event = runmgr->GetCurrentEvent();
		int evId = event->GetEventID();
		G4ThreeVector position = track->GetPosition();
		G4double time = track->GetGlobalTime();
		G4ThreeVector momentum = track->GetMomentum();

		// get B and E fields
		G4ThreeVector B, E;
		if(fields) {
			G4double point[4], field[6];
			point[0] = position[0];
			point[1] = position[1];
			point[2] = position[2];
			point[3] = time;
			BLGlobalField::getObject()->GetFieldValue(point,field);
			B = G4ThreeVector(field[0],field[1],field[2]);
			E = G4ThreeVector(field[3],field[4],field[5]);
		}

		// transform to desired coordinates, if available
		BLCoordinates *c = (BLCoordinates *)track->GetUserInformation();
		if(c && c->isValid()) {
			c->getCoords(coordType,position);
			momentum = c->getRotation() * momentum;
			B = c->getRotation() * B;
			E = c->getRotation() * E;
		}
		setVariable((G4String("x")+suffix).c_str(), position[0]/mm);
		setVariable((G4String("y")+suffix).c_str(), position[1]/mm);
		setVariable((G4String("z")+suffix).c_str(), position[2]/mm);
		setVariable((G4String("Px")+suffix).c_str(), momentum[0]/MeV);
		setVariable((G4String("Py")+suffix).c_str(), momentum[1]/MeV);
		setVariable((G4String("Pz")+suffix).c_str(), momentum[2]/MeV);
		setVariable((G4String("t")+suffix).c_str(), time/ns);
		setVariable((G4String("PDGid")+suffix).c_str(), track->GetDefinition()->GetPDGEncoding());
		setVariable((G4String("EventID")+suffix).c_str(), evId);
		setVariable((G4String("TrackID")+suffix).c_str(), 
				BLManager::getObject()->
				    getExternalTrackID(track));
		setVariable((G4String("ParentID")+suffix).c_str(),
				BLManager::getObject()->
				    getExternalParentID(track));
		setVariable((G4String("wt")+suffix).c_str(),track->GetWeight());
		if(fields) {
			setVariable((G4String("Bx")+suffix).c_str(),B[0]/tesla);
			setVariable((G4String("By")+suffix).c_str(),B[1]/tesla);
			setVariable((G4String("Bz")+suffix).c_str(),B[2]/tesla);
			setVariable((G4String("Ex")+suffix).c_str(),
							E[0]/(megavolt/meter));
			setVariable((G4String("Ey")+suffix).c_str(),
							E[1]/(megavolt/meter));
			setVariable((G4String("Ez")+suffix).c_str(),
							E[2]/(megavolt/meter));
		}
	}
#endif // NO_GEANT4

	/// isidchar() returns true if the char is valid for an identifier.
	bool isidchar(char c) { return isalnum(c) || c == '_'; }

	/// isopchar() returns true if the character is an operator (no parens)
	bool isopchar(char c) { return c=='+' || c=='-' || c=='*' || c=='/' ||
					c=='^' || c=='<' || c=='>' || c=='!' ||
					c=='|' || c=='&' || c=='=' || c=='%'; }

	/// isumchar() returns true if the char should be inside the parens
	/// generated to surround the unary minus (parens matched in evaluate())
	bool isumchar(char c) { return isidchar(c) || c=='.' ||
					c=='^' || c=='(' || c==')'; }
};

#endif // BLEVALUATOR_HH
