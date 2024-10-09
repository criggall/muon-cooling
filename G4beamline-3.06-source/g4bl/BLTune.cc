//	BLTune.cc

#include "BLTune.hh"
#include "BLCommand.hh"

BLEvaluator BLTune::eval;
std::map<G4double*,G4String> BLTune::tuneExpr;
std::map<G4double*,G4double> BLTune::tuneUnits;
std::set<G4String> BLTune::names;

void BLTune::defineTunableArg(G4double& arg, G4double units, G4String expr)
{
	tuneUnits[&arg] = units;
	tuneExpr[&arg] = expr;
	if(!eval.findFunction("sin",1)) eval.setStdMath();
	arg = eval.evaluate(expr) * units;
	if(eval.status() != HepTool::Evaluator::OK) {
		BLCommand::printError("Invalid expression '%s'",expr.c_str());
	}
}

void BLTune::copyTunableArg(G4double *newArg, const G4double *oldArg)
{
	if(tuneExpr.count((G4double*)oldArg) == 0)
		defineTunableArg(*newArg,1,"0");
	else
		defineTunableArg(*newArg,tuneUnits[(G4double*)oldArg],
						tuneExpr[(G4double*)oldArg]);
	*newArg = *oldArg;
}

bool BLTune::isSet(G4String name)
{
	return eval.findVariable(name);
}

void BLTune::set(G4String& name, G4double value)
{
	names.insert(name);
	eval.setVariable(name,value);
	update();
}

void BLTune::unset(G4String& name)
{
	eval.removeVariable(name);
}

void BLTune::update()
{
	if(!eval.findFunction("sin",1)) eval.setStdMath();

	std::map<G4double*,G4String>::iterator i;
	for(i=tuneExpr.begin(); i!=tuneExpr.end(); ++i) {
		G4double *p = i->first;
		G4String expr = i->second;
		G4double units = tuneUnits[p];
		*p = eval.evaluate(expr) * units;
		if(eval.status() != HepTool::Evaluator::OK) {
			BLCommand::printError("Invalid expression '%s'",
							expr.c_str());
		}
	} 
}
