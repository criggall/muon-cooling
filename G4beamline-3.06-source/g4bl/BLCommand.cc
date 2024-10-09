//	BLCommand.cc
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

#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <cmath>	// std::isnan()

#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4UImanager.hh"

#include "BLCommand.hh"
#include "G4NistManager.hh"
#include "BLParam.hh"
#include "BLTune.hh"
#include "BLEvaluator.hh"
#include "BLNTuple.hh"
#include "BLSignal.hh"
#include "BLMPI.hh"
#include "fnmatch.h"

#ifdef WIN32
#define snprintf _snprintf
#endif
#ifdef __CYGWIN__
#include "mysnprintf.hh"
#endif


// mapCommand cannot be a simple object, as it is required in registerCommand(),
// which is called from static initializers, and we cannot guarantee the order.
// registerCommand() creates the map the first time it is called.
std::map<G4String,BLCommand*> *BLCommand::mapCommand;

int BLCommand::errors = 0;
std::istream *BLCommand::in = 0;

BLCommand::~BLCommand()
{
}

void BLCommand::registerCommand(BLCmdType _type)
{
	type = _type;

	if(!mapCommand)
		mapCommand = new std::map<G4String,BLCommand*>;
	if((*mapCommand).count(commandName()) > 0 &&
	   (*mapCommand)[commandName()] != this)
		printError("Attempt to re-define the command '%s'",
				commandName().c_str());
	else
		(*mapCommand)[commandName()] = this;
}

void BLCommand::printError(const char *fmt, ...)
{
	va_list ap;
	va_start(ap,fmt);
	char tmp[1024];
	vsprintf(tmp,fmt,ap);
	G4Exception("","printError",JustWarning,tmp);
	++errors;
}

G4String BLCommand::wrapWords(G4String text, G4String indent1,
				G4String indent, G4String::size_type width)
{
	G4String retval(""), space(" "), line(""), word("");
	G4String::size_type pos=0, next=0;
	bool eol=false;

	while(pos < text.size()) {
		if(eol && isspace(text[pos])) {
			next = text.find_first_of("\n",pos);
			if(next == text.npos) next = text.size();
			if(line.size() > 0) {retval += line; retval += "\n";}
			line = indent + text.substr(pos,next-pos);
			retval += line;
			retval += "\n";
			line = "";
			pos = next + 1;
			// eol = true;
			continue;
		}
		next = text.find_first_of(" \t\r\n",pos);
		if(next == text.npos) {
			word = text.substr(pos);
			next = text.size();
		} else {
			word = text.substr(pos,next-pos);
			eol = (text[next] == '\n');
		}
		if(word.size() > 0) {
			if(line.size() == 0) {
				if(retval.size() == 0)
					line = indent1;
				else
					line = indent;
				line += word;
			} else if(line.size()+word.size()+1 < width) {
				line += space;
				line += word;
			} else {
				retval += line;
				retval += "\n";
				line = indent;
				line += word;
			} 
		}
		pos = next + 1;
	} 
	if(line.size() > 0) {retval += line; retval += "\n";}

	return retval;
}

void BLCommand::argString(G4String& var, G4String name, G4String _description,
		bool permitChange)
{
	switch(argMode) {
	case CHANGE:
		if(!permitChange) {
			if(argName == name)
				printError("%s: Argument '%s' is not permitted to change\n",
					commandName().c_str(),name.c_str());
			break;
		}
		// flow into PROCESS
	case PROCESS:	
		if(argFound || argName != name) return;
		argFound = true;
		var = argValue;
		break;
	case COUNT:
		++nArgs;
		if(!permitChange) ++nFixed;
		break;
	case HELP:	
		printArgDesc(name,_description+(permitChange?"":" #"));
		break;
	case PRINT:
		printArgString += name;
		printArgString += "=";
		printArgString += var;
		printArgString += "\n";
		break;
	}
}

void BLCommand::argDouble(G4double& var, G4String name, G4String _description,
			G4double units, G4String fmt, bool permitChange)
{
	if(fmt == "") fmt = DefaultDoubleFmt;

	switch(argMode) {
	case CHANGE:
		if(!permitChange) {
			if(argName == name)
				printError("%s: Argument '%s' is not permitted to change\n",
					commandName().c_str(),name.c_str());
			break;
		}
		// flow into PROCESS
	case PROCESS:	
		if(argFound || argName != name) return;
		argFound = true;
		{ BLEvaluator e;
		  G4double v = e.evaluate(argValue);
		  if(e.status() != HepTool::Evaluator::OK) {
			  printError("Invalid value '%s' for argument '%s'"
				" to command '%s'", argValue.c_str(),
				argName.c_str(), commandName().c_str());
			  return;
		  }
		  var = v * units;
		}
/**** 
		{ const char *p=argValue.c_str();
		  char *q;
		  G4double v = strtod(p,&q);
		  if(p == q || *q != '\0') {
			  printError("Invalid value '%s' for argument '%s'"
				" to command '%s'", argValue.c_str(),
				argName.c_str(), commandName().c_str());
		  	return;
		  }
		  var = v * units;
		}
****/
		break;
	case COUNT:
		++nArgs;
		if(!permitChange) ++nFixed;
		break;
	case HELP:	
		printArgDesc(name,_description+(permitChange?"":" #"));
		break;
	case PRINT:
		{ char tmp[65];
		  snprintf(tmp,sizeof(tmp),fmt.c_str(),var/units);
		  printArgString += name;
		  printArgString += "=";
		  printArgString += tmp;
		  printArgString += "\n";
		}
		break;
	}
}

void BLCommand::argTunable(G4double& var, G4String name, G4String _description,
			G4double units, G4String fmt)
{
	if(fmt == "") fmt = DefaultDoubleFmt;

	switch(argMode) {
	case CHANGE:
	case PROCESS:	
		if(argFound || argName != name) return;
		argFound = true;
		BLTune::defineTunableArg(var,units,argValue);
		break;
	case HELP:	
		printArgDesc(name,_description+" @");
		break;
	case COUNT:
		++nArgs;
		++nTunable;
		break;
	case PRINT:
		{ char tmp[65];
		  snprintf(tmp,sizeof(tmp),fmt.c_str(),var/units);
		  printArgString += name;
		  printArgString += "=";
		  printArgString += tmp;
		  printArgString += "\n";
		}
		break;
	}
}

void BLCommand::argInt(G4int& var, G4String name, G4String _description,
		bool permitChange)
{
	switch(argMode) {
	case CHANGE:
		if(!permitChange) {
			if(argName == name)
				printError("%s: Argument '%s' is not permitted to change\n",
					commandName().c_str(),name.c_str());
			break;
		}
		// flow into PROCESS
	case PROCESS:	
		if(argFound || argName != name) return;
		argFound = true;
		{ BLEvaluator e;
		  G4double v = e.evaluate(argValue);
		  var = (G4int)(v+0.5);
		  if(e.status() != HepTool::Evaluator::OK ||
							fabs(var-v) > 0.001) {
			  printError("Invalid value '%s' for argument '%s'"
				" to command '%s'", argValue.c_str(),
				argName.c_str(), commandName().c_str());
			  return;
		  }
		}
/***
		{ const char *p=argValue.c_str();
		  char *q;
		  G4int v = strtol(p,&q,0);
		  if(p == q || *q != '\0') {
			  printError("Invalid value '%s' for argument '%s'"
				" to command '%s'", argValue.c_str(),
				argName.c_str(), commandName().c_str());
		  	return;
		  }
		  var = v;
		}
***/
		break;
	case COUNT:
		++nArgs;
		if(!permitChange) ++nFixed;
		break;
	case HELP:	
		printArgDesc(name,_description+(permitChange?"":" #"));
		break;
	case PRINT:
		{ char tmp[65];
		  snprintf(tmp,sizeof(tmp),"%d",var);
		  printArgString += name;
		  printArgString += "=";
		  printArgString += tmp;
		  printArgString += "\n";
		}
		break;
	}
}

void BLCommand::printArgDesc(G4String name, G4String _description)
{
	G4String indent1(IndentDesc);
	indent1 += name;
	do { indent1 += " "; } while(indent1.size() < IndentArg.size());
	printf("%s",wrapWords(_description,indent1,IndentArg).c_str());
}

int BLCommand::handleNamedArgs(BLArgumentMap& args)
{
	int retval = 0;
	argMode = PROCESS;

	BLArgumentMap::iterator i;
	for(i=args.begin(); i!=args.end(); ++i) {
		argName = i->first;
		argValue = i->second;
		argFound = false;
		defineNamedArgs();
		if(!argFound) { 
			printError("Invalid argument '%s' to command '%s'",
				argName.c_str(),commandName().c_str());
			retval = -1;
		}
	}

	if(args.size() > 0) argChanged();

	return retval;
}

void BLCommand::defineNamedArgs()
{
	if(argMode == HELP) printf("%s(none)\n",IndentDesc.c_str());
}

void BLCommand::help(bool detailed)
{
	G4String str(commandName());
	do { str += " "; } while(str.size() < IndentDesc.size());
	printf("%s%s\n",str.c_str(),synopsis.c_str());
	if(detailed) {
		if(description.rfind(": ") == description.size()-2)
			description += BLNTuple::getFormatList(); 
		printf("\n%s",wrapWords(description,IndentDesc,IndentDesc).c_str());
		nArgs = nFixed = nTunable = 0;
		argMode = COUNT;
		defineNamedArgs();
		if(nArgs > 0) {
			printf("\n%sNamed Arguments", IndentDesc.c_str());
			if(nFixed > 0) 
				printf(" (#=cannot be changed in place cmd)");
			if(nTunable > 0) printf(" (@=Tunable)");
			printf(":\n");
			argMode = HELP;
			defineNamedArgs();
		}
	}
}

void BLCommand::printArgs(G4String indent1) 
{
	printArgString = "";
	argMode = PRINT;
	defineNamedArgs();
	printf("%s",wrapWords(printArgString,indent1,IndentArg).c_str());
}

void BLCommand::print(G4String name)
{
	G4String indent1(commandName());
	do { indent1 += " "; } while(indent1.size() < IndentDesc.size());
	indent1 += name;
	do { indent1 += " "; } while(indent1.size() < IndentArg.size());
	printArgs(indent1);
}

void BLCommand::print(G4String name, BLArgumentMap& namedArgs)
{
	G4String indent1(commandName());
	do { indent1 += " "; } while(indent1.size() < IndentDesc.size());
	indent1 += name;
	do { indent1 += " "; } while(indent1.size() < IndentArg.size());
	printArgString = "";
	BLArgumentMap::iterator i;
	for(i=namedArgs.begin(); i!=namedArgs.end(); ++i) {
		printArgString += i->first;
		printArgString += "=";
		printArgString += i->second;
		printArgString += "\n";
	}
	printf("%s",wrapWords(printArgString,indent1,IndentArg).c_str());
}

int BLCommand::doCommand(G4String& line)
{
	BLMPI::scan();

	BLArgumentVector argv;
	BLArgumentMap namedArgs;
	G4String::size_type place = 0;
	TokenType type;

	if(line.size() < 1) return 0;
	char c = line.c_str()[0];
	if(c == '#') return 0;
	if(c == '*') {
		printf("%s\n",line.c_str());
		return 0;
	}
	if(c == '!') {
		printf("%s\n",line.c_str());
		fflush(stdout);
		system(line.c_str()+1);
		return 0;
	}
	if(c == '/') {
		printf("%s\n",line.c_str());
		G4UImanager* UI = G4UImanager::GetUIpointer();
		UI->ApplyCommand(line.c_str());
		return 0;
	}

	G4String cmd = nextToken(line,place,type);
	switch(type) {
	case NONE:
		return 0;	// line with only whitespace
	case ARGVALUE:
		break;
	case ARGNAME:
		printError("Syntax Error 1");
		return -1;
	}

	if(parseArgs(line.substr(place),argv,namedArgs) < 0)
		return -1;

	int retval;
	if((*mapCommand).count(cmd) > 0) {
		retval = (*mapCommand)[cmd]->command(argv,namedArgs);
	} else {
		retval = -1;
		printError("Unknown command '%s'",cmd.c_str());
	}

	return retval;
}

int BLCommand::parseArgs(const G4String &line, BLArgumentVector &argv,
						BLArgumentMap &namedArgs)
{
	G4String::size_type place = 0;
	while(place < line.size()) {
		TokenType type;
		G4String arg = nextToken(line,place,type), val;
		switch(type) {
		case NONE:
			break;
		case ARGNAME:
			val = nextToken(line,place,type);
			if(type != ARGVALUE && !(type == NONE && 
							place >= line.size())) {
				printError("Syntax Error 2");
				return -1;
			}
			namedArgs[arg] = Param.expand(val);
			break;
		case ARGVALUE:
			argv.push_back(Param.expand(arg));
			break;
		}
	}
	return 0;
}

G4String BLCommand::nextToken(const G4String& line, G4String::size_type& place, 
							TokenType& type)
{
	// (add +- for particlecolor pi+=1,1,1)
	static const char namechars[] = ",+-ABCDEFGHIJKLMNOPQRSTUVWXYZ_"
				       "abcdefghijklmnopqrstuvwxyz0123456789";
	G4String::size_type i;

	// check if previous token was ARGNAME
	if(line[place] == '=') {
		++place;
		goto value;
	}

	// skip initial whitespace
	while(place < line.size() && isspace(line[place])) ++place;

	// check for End of line or comment-to-end-of-line
	if(place >= line.size() || line[place] == '#') {
		type = NONE;
		place = line.size();
		return "";
	}

	// check for ARGNAME
	if(isalnum(line[place]) || line[place] == '_') {
		i = line.find_first_not_of(namechars,place);
		if(i > place && i < line.size() && line[i] == '=' &&
							line[i+1] != '=') {
			G4String retval = line.substr(place,i-place);
			place = i;
			type = ARGNAME;
			return retval;
		}
	}
value:
	if(line[place] == '"') {
		++place;
		i = line.find('"',place);
		if(i <line.size()) {
			G4String retval = line.substr(place,i-place);
			place = i + 1;
			type = ARGVALUE;
			return retval;
		}
	} else if(line[place] == '\'') {
		++place;
		i = line.find('\'',place);
		if(i <line.size()) {
			G4String retval = line.substr(place,i-place);
			place = i + 1;
			type = ARGVALUE;
			return retval;
		}
	}

	if(place >= line.size()) {
		type = NONE;
		return "";
	}

	// find next whitespace
	G4String::size_type start = place;
	while(place < line.size() && !isspace(line[place])) ++place;

	type = ARGVALUE;

	return line.substr(start,place-start);
}

int BLCommand::readFile(G4String filename)
{
	std::istream *save = in;

	std::ifstream fin;
	if(filename == "-") {
		in = &std::cin;
	} else {
#ifdef WIN32
		fin.open(filename.c_str(),fin.in|fin.binary);
#else
		fin.open(filename.c_str());
#endif
		in = &fin;
	}
	if(in == 0 || !in->good()) {
		printError("Cannot open file '%s'",filename.c_str());
		return -1;
	}

	int errorCount = getErrorCount();

	while(in->good()) {
		G4String *line = getNextCommand();
		if(!line || line->find("exit") == 0) break;
		doCommand(*line);
	}

	fin.close();
	if(in == &std::cin) std::cout << std::endl;

	in = save;
	return getErrorCount() - errorCount;
}

G4String *BLCommand::getNextCommand()
{
	static G4String line;

	if(in == &std::cin) {
		BLSignal::setSignalReceived(true);
		std::cout << "cmd: ";
	}

	char tmp[1024+1];

	do {
		if(!in->good()) return 0;
		in->getline(tmp,sizeof(tmp));
		// trim trailing whitespace (handle "...\\ \r" from Windows)
		for(int i=strlen(tmp)-1; i>=0&&isspace(tmp[i]); --i)
			tmp[i] = '\0';
		line = tmp;
	} while(line.size() == 0);

	// handle continuation lines (ending with \)
	while(line[line.size()-1] == '\\' && in->good()) {
		line[line.size()-1] = ' ';
		in->getline(tmp,sizeof(tmp));
		// trim trailing whitespace (handle "...\\ \r" from Windows)
		for(int i=strlen(tmp)-1; i>=0&&isspace(tmp[i]); --i)
			tmp[i] = '\0';
		line += tmp;
	}

	if(in == &std::cin)
		BLSignal::setSignalReceived(false);

	// trim whitespace
	G4String::size_type i = line.find_first_not_of(" \t\r\n\v");
	if(i == line.npos) i = line.size();
	line.erase(0,i);
	while((i=line.find_last_of(" \t\r\n\v")) != line.npos) {
		if(i != line.size()-1) break;
		line.erase(i);
	}

	return &line;
}

BLCommandPos BLCommand::getPos()
{
	BLCommandPos pos(in);
	return pos;
}

void BLCommand::setPos(BLCommandPos &pos)
{
	if(in != pos.in) {
		printError("Invalid BLCommand setPos()");
		return;
	}

	in->seekg(pos.pos);
}

G4Material *BLCommand::getMaterial(G4String materialName, bool ignoreError)
{
	const G4String G4("G4_");
	G4NistManager *m = G4NistManager::Instance();

	// first, check if it has already been constructed
	G4Material *mat = G4Material::GetMaterial(materialName,false);
	if(mat) return mat;
	mat = G4Material::GetMaterial(G4+materialName,false);
	if(mat) return mat;

	// Check for alias
	// NOTE: these materials are also listed in BLCMDmaterial::help()
	if(materialName == "Air") materialName = "AIR";
	if(materialName == "H2O") materialName = "WATER";
	if(materialName == "LH2") materialName = "lH2";
	if(materialName == "Vacuum") {
		// preserve the name "Vacuum" 
		G4Material *g = m->FindOrBuildMaterial("G4_Galactic");
		G4Material *v = new G4Material(materialName,
			g->GetDensity(),1,g->GetState(),g->GetTemperature(),
			g->GetPressure());
		v->AddMaterial(g,1.0);
		return v;
	}
	if(materialName == "Scintillator") {
		// preserve the name "Scintillator" 
		G4Material *g = m->FindOrBuildMaterial("G4_POLYSTYRENE");
		G4Material *v = new G4Material(materialName,
			g->GetDensity(),1,g->GetState(),g->GetTemperature(),
			g->GetPressure());
		v->AddMaterial(g,1.0);
		return v;
	}
	if(materialName == "Stainless304") {
		// ignore optional components < 1%
		G4Material *v = new G4Material(materialName,8.01*g/cm3,4,
			kStateSolid);
		G4Material *Cr = m->FindOrBuildMaterial("G4_Cr");
		G4Material *Ni = m->FindOrBuildMaterial("G4_Ni");
		G4Material *Mn = m->FindOrBuildMaterial("G4_Mn");
		G4Material *Fe = m->FindOrBuildMaterial("G4_Fe");
		v->AddMaterial(Cr,0.19);
		v->AddMaterial(Ni,0.10);
		v->AddMaterial(Mn,0.01);
		v->AddMaterial(Fe,0.70);
		return v;
	}
	if(materialName == "Stainless316") {
		// ignore optional components < 1%
		G4Material *v = new G4Material(materialName,8.01*g/cm3,5,
			kStateSolid);
		G4Material *Cr = m->FindOrBuildMaterial("G4_Cr");
		G4Material *Ni = m->FindOrBuildMaterial("G4_Ni");
		G4Material *Mn = m->FindOrBuildMaterial("G4_Mn");
		G4Material *Mo = m->FindOrBuildMaterial("G4_Mo");
		G4Material *Fe = m->FindOrBuildMaterial("G4_Fe");
		v->AddMaterial(Cr,0.18);
		v->AddMaterial(Ni,0.12);
		v->AddMaterial(Mn,0.01);
		v->AddMaterial(Mo,0.025);
		v->AddMaterial(Fe,0.665);
		return v;
	}
	if(materialName == "lHe" || materialName == "LHe") {
		G4Material *v = new G4Material(materialName,0.1249*g/cm3,1,
			kStateLiquid);
		G4Material *He = m->FindOrBuildMaterial("G4_He");
		v->AddMaterial(He,1.0);
		return v;
	}

	// all materials in the database begin "G4_".
	if(materialName.find(G4) == 0)
		mat = m->FindOrBuildMaterial(materialName);
	else
		mat = m->FindOrBuildMaterial(G4+materialName);
	if(mat) return mat;

	if(ignoreError) {
		return 0;
	}

	G4Exception("getMaterial","Material not found",FatalException,
							materialName.c_str());
	return 0;
}

G4VisAttributes *BLCommand::getVisAttrib(G4String color)
{
	G4VisAttributes *p=0;

	std::vector<G4double> v = getList(color,',');
	if(v.size() == 3) {
		p = new G4VisAttributes(true,G4Color(v[0],v[1],v[2]));
	} else if(v.size() == 4) {
		p = new G4VisAttributes(true,G4Color(v[0],v[1],v[2],v[3]));
	} else {
		p = new G4VisAttributes(false); // invisible
		if(color.size() > 0 && toupper(color.c_str()[0]) != 'I') {
			printError("Invalid color '%s'",color.c_str());
		}
	}

	p->SetDaughtersInvisible(false);

	return p;
}

G4RotationMatrix *BLCommand::stringToRotationMatrix(G4String rotation)
{
	// We apply successive rotations OF THE OBJECT around the FIXED
	// axes of the parent's local coordinates; rotations are applied
	// left-to-right (rotation="r1,r2,r3" => apply r1 then r2 then r3).
	G4RotationMatrix *rot = new G4RotationMatrix();
	std::vector<G4String> v = splitString(rotation,",",true);
	if(v.size() == 4 && toupper(v[0].c_str()[0]) == 'A') {
		double x = getDouble(v[0].substr(1));
		double y = getDouble(v[1]);
		double z = getDouble(v[2]);
		double angle = getDouble(v[3]) * deg;
		G4ThreeVector axis(x,y,z);
		*rot = G4RotationMatrix(axis.unit(),angle) * *rot;
	} else {
	    for(unsigned i=0; i<v.size(); ++i) {
		G4double angle=getDouble(v[i].substr(1)) * deg;
		if(std::isnan(angle)) {
			printError("Invalid rotation specification '%s'",
						rotation.c_str());
			return 0;
		}
		G4RotationMatrix thisRotation;
		switch(v[i].c_str()[0]) {
		case 'X': case 'x':
			thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
			break;
		case 'Y': case 'y':
			thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
			break;
		case 'Z': case 'z':
			thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
			break;
		default:
			printError("Invalid rotation specification '%s'",
						rotation.c_str());
			return 0;
		}
		*rot = thisRotation * *rot;
	    }
	}

	return rot;
}

void BLCommand::dumpRotation(const G4RotationMatrix *rot, const char *str)
{
	if(!rot) {
		printf("%s: no rotation\n",str);
		return;
	}
	G4ThreeVector x(1.0,0.0,0.0), z(0.0,0.0,1.0);
	G4ThreeVector rx = *rot * x;
	G4ThreeVector rz = *rot * z;
	printf("%s: %.3f,%.3f,%.3f / %.3f,%.3f,%.3f / %.3f,%.3f,%.3f\n",
		str,rot->xx(),rot->xy(),rot->xz(),rot->yx(),rot->yy(),rot->yz(),
		rot->zx(),rot->zy(),rot->zz());
	printf("        rot*x=%.3f,%.3f,%.3f   rot*z = %.3f,%.3f,%.3f\n",
		rx[0],rx[1],rx[2],rz[0],rz[1],rz[2]);
}

std::vector<G4String> BLCommand::splitString(const G4String &s, 
					const G4String &delim, bool trim)
{
	G4String tmp(s), p;
	std::vector<G4String> v;
	while(tmp.size() > 0) {
		G4String::size_type i=tmp.find_first_of(delim);
		if(i > tmp.size()) i = tmp.size();
		p = tmp.substr(0,i);
		if(trim) {
			while(p.size() > 0 && isspace(p.c_str()[0]))
				p.erase(0,1);
			while(p.size() > 0 && isspace(p.c_str()[p.size()-1]))
				p.erase(p.size()-1,1);
		}
		v.push_back(p);
		tmp.erase(0,i+1);
	}

	return v;
}

std::vector<G4double> BLCommand::getList(const G4String &s,
							const G4String &delim)
{
	std::vector<G4String> vect_s = splitString(s,delim,true);
	std::vector<G4double> retval;

	for(unsigned i=0; i<vect_s.size(); ++i) {
		G4double val=getDouble(vect_s[i]);
		if(std::isnan(val)) {
			retval.clear();
			break;
		}
		retval.push_back(val);
	}

	return retval;
}

bool BLCommand::matchList(const G4String s, const G4String patternList)
{
	std::vector<G4String> v = splitString(patternList,",",false);
	for(unsigned i=0; i<v.size(); ++i) {
		if(fnmatch(v[i].c_str(),s.c_str(),0) == 0)
			return true;
	}

	return false;
}

double BLCommand::undefined()
{
	return std::strtod("nan()",0);
}

bool BLCommand::isUndefined(double v)
{
	return std::isnan(v);
}
