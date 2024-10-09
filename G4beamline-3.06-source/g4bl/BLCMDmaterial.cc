//	BLCMDmaterial.cc
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

#include "globals.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ParticleTable.hh"

#include "BLManager.hh"
#include "BLParam.hh"
#include "BLCommand.hh"
#include "BLEvaluator.hh"
#include "G4NistManager.hh"

#ifdef _MSC_VER
#define snprintf _snprintf
#endif
#ifdef __CYGWIN__
#include "mysnprintf.hh"
#endif

const G4double UNDEFINED = -1.0e30;

/**	class BLCMDmaterial implements the material command, which
 *	defines materials.
 *
 **/
class BLCMDmaterial : public BLCommand
{
	G4double a;
	G4double z;
	G4double density;
	G4double pressure;
	G4double temperature;
	G4String state;
	G4String keep;
	G4String kill;
	G4String require;
	G4Material *material;
	bool complete_description;
public:
	/// Constructor.
	BLCMDmaterial();

	/// Destructor.
	~BLCMDmaterial() { }

	/// commandName() returns "material"
	G4String commandName()  { return "material"; }

	/// command() executes the command associated with this element.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command. 
	void defineNamedArgs();

	/// help() prints the help.
	void help(bool detailed);
};
BLCMDmaterial defaultMaterial;

class MaterialFilter : public BLManager::SteppingAction, public BLCallback
{
	std::set<G4ParticleDefinition*> keepSet;
	std::set<G4ParticleDefinition*> killSet;
	BLEvaluator *eval;
	G4String require;
	long nKilled;
	BLManager *manager;
	G4Material *material;
	G4ParticleTable *table;
public:
	MaterialFilter(G4String keep, G4String kill, G4String _require,
						G4Material *_material);
	void UserSteppingAction(const G4Step *step);
	void callback(int type);
};

BLCMDmaterial::BLCMDmaterial() : BLCommand()
{
	registerCommand(BLCMDTYPE_AUX);
	setSynopsis("construct a new material.");
	setDescription("This is an interface to G4Material.\n"
		"This command is rarely required, because elements and most "
		"common materials are available via the NIST database. "
		"Any material available from the NIST database can simply "
		"be used -- if it is unknown then it "
		"will be automatically defined from the database. "
		"Uncommon materials or nonstandard densities must be defined "
		"with this command.\n\n"
		"The first argument to this command is the material name, "
		"which is always required; density is also required. "
		"The command to define an element (e.g. with non-standard "
		"density) is:\n"
		"    material H2 Z=1 A=1.01 density=0.000090\n"
		"A mixture or compound is a combination of known materials "
		"and/or elements; the command is:\n"
		"    material water H,0.1119 O,0.8881 density=1.0\n"
		"The numbers following the element names are their mass "
		"fractions (note that WATER is available from the NIST db). "
		"Either type of command can optionally have: pressure, "
		"temperature, state.\n\n"
		"With no arguments, this command prints the current material "
		"table. Note that 'G4_' is prepended to the names of most "
		"materials that are obtained from the NIST database; 'G4_Al' "
		"and 'Al' refer to the same material (unless one was "
		"previously defined using this command).\n\n"
		"The following three arguments permit track filtering for all "
		"volumes made of this material:\n"
		"  keep    A comma-separated list of particle names to keep.\n"
		"  kill    A comma-separated list of particle names to kill.\n"
		"  require An expression that must evaluate nonzero or the "
		"track\n"
		"          is killed.\n"
		"The require expression uses global coordinates and can use "
		"the followng track variables:\n"
		"  x,y,z,Px,Py,Pz,t,PDGid,EventID,TrackID,ParentID,wt\n\n"
		"The following materials are known from the NIST database, "
		"and will be automatically created on first use:\n\n");
	// cannot call GetNistMaterialNames() here, so defer to the help()
	// function, based on complete_description.

	// default values
	a = z = density = UNDEFINED;
	pressure = STP_Pressure;
	temperature = STP_Temperature;
	state = "";
	keep = "";
	kill = "";
	require = "";
	material = 0;
	complete_description = false;
}

int BLCMDmaterial::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 1) {
		G4cout << *G4Material::GetMaterialTable() << G4endl;
		return 0;
	}

	// default values
	a = z = density = UNDEFINED;
	pressure = STP_Pressure;
	temperature = STP_Temperature;
	state = "";
	keep = "";
	kill = "";
	require = "";
	material = 0;

	G4String name = argv[0];

	int retval = handleNamedArgs(namedArgs);

	G4State s = kStateUndefined;
	switch(state.c_str()[0]) {
	case 'g': case 'G':	s = kStateGas;		break;
	case 'l': case 'L':	s = kStateLiquid;	break;
	case 's': case 'S':	s = kStateSolid;	break;
	}

	if(density == UNDEFINED) {
		printError("material: density is required.");
		return -1;
	}

	if(argv.size() > 1) {		// mixture of other materials
		int nComponents = argv.size() - 1;
		material = new G4Material(name,density,nComponents,s,
						temperature,pressure);
		int i;
		G4double sum = 0.0;
		G4String mix;
		for(i=1; i<=nComponents; ++i) {
			G4String::size_type j = argv[i].find(',');
			if(j == argv[i].npos) break;
			G4String p = argv[i].substr(0,j);
			G4Material *m = getMaterial(p);
			if(!m) break;
			char *q = 0;
			G4double f = strtod(argv[i].substr(j+1).c_str(),&q);
			if(!q || *q != '\0') break;
			material->AddMaterial(m,f);
			sum += f;
			char tmp[64];
			snprintf(tmp,sizeof(tmp),"%.2f*%s",f,p.c_str());
			if(i > 1) mix += "+";
			mix += tmp;
		}
		if(i != nComponents+1 || sum < 0.99 || sum > 1.01) {
			printError("material %-6s invalid mixture/compound",
				name.c_str());
			return -1;
		}
		printf("material %-6s Mixture: %s\n",name.c_str(),mix.c_str());
		printf("                density=%.3f temperature=%.0f pressure=%.1f\n",
			density/(gram/cm3),temperature,pressure/atmosphere);
	} else {			// standalone material
		if(a == UNDEFINED || z == UNDEFINED) {
			printError("material: need A and/or Z");
			return -1;
		}
		material = new G4Material(name,density,1,s,
						temperature,pressure);
		G4Element *elem = new G4Element(name,name,1);
		elem->AddIsotope(new G4Isotope(name,(int)(z+0.5),
						(int)(a/gram*mole+0.5)),1.0);
		material->AddElement(elem,1.0);
		printf("material %-6s Z=%.2f A=%.2f\n",name.c_str(),z,
							a/(gram/mole));
		printf("                density=%.3f temperature=%.0f pressure=%.1f\n",
			density/(gram/cm3),temperature,pressure/atmosphere);
	}

	if(keep != "" || kill != "" || require != "")
		new MaterialFilter(keep,kill,require,material);

	return retval;
}

void BLCMDmaterial::defineNamedArgs()
{
	argDouble(a,"a","Effective Atomic Weight of the material (g/mole)",gram/mole);
	argDouble(z,"z","Effective Atomic Number of the material");
	argDouble(density,"density","Density of the material (gm/cm^3)",gram/cm3);
	argDouble(pressure,"pressure","Pressure of the material (Atm)",atmosphere);
	argDouble(temperature,"temperature","Temperature of the material (K)");
	argString(state,"state","State of the material (g,l, or s)");
	argDouble(a,"A","Synonym for a (g/mole)",gram/mole);
	argDouble(z,"Z","Synonym for z");
	argString(keep,"keep","A comma-separated list of particle names to keep"
				" (all others are killed; ignored if empty).");
	argString(kill,"kill","A comma-separated list of particle names to kill.");
	argString(require,"require","An expression that must evaluate nonzero "
			"or the track is killed.");
	argDouble(a,"aaa","Effective Atomic Weight of the material (g/mole)",gram/mole);

}

void BLCMDmaterial::help(bool detailed)
{
	if(!complete_description) {
		complete_description = true;
		G4NistManager *m = G4NistManager::Instance();
		const std::vector<G4String> &v=m->GetNistMaterialNames();
		for(unsigned i=0; i<v.size(); ++i) {
			const char *p = v[i].c_str();
			if(strncmp(p,"G4_",3) == 0) p += 3; // omit the initial "G4_"
			description += p;
			description += " ";
		}
		description += "Stainless304 Stainless316 lHe ";
		description += "\n\nAliases: LHe=lHe Air=AIR, H2O=WATER, "
				"Vacuum=Galactic, LH2=lH2, "
				"Scintillator=POLYSTYRENE\n";
	}

	BLCommand::help(detailed);
}

MaterialFilter::MaterialFilter(G4String keep, G4String kill, G4String _require,
			G4Material *_material) :  BLManager::SteppingAction(),
					BLCallback(), keepSet(), killSet()
{
	eval = 0;
	require = _require;
	nKilled = 0;
	manager = BLManager::getObject();
	material = _material;
	table = G4ParticleTable::GetParticleTable();

	if(keep.size() != 0) {
		printf("                keep=%s\n",keep.c_str());
		std::vector<G4String> v=BLCommand::splitString(keep,',',true);
		for(unsigned i=0; i<v.size(); ++i) {
			G4ParticleDefinition *p = table->FindParticle(v[i]);
			if(!p) {
				BLCommand::printError("material: particle '%s' not found",
					v[i].c_str());
			}
			keepSet.insert(p);
		}
	}

	if(kill.size() != 0) {
		printf("                kill=%s\n",kill.c_str());
		std::vector<G4String> v=BLCommand::splitString(kill,',',true);
		for(unsigned i=0; i<v.size(); ++i) {
			G4ParticleDefinition *p = table->FindParticle(v[i]);
			if(!p) {
				BLCommand::printError("material: particle '%s' not found",
					v[i].c_str());
			}
			killSet.insert(p);
		}
	}

	if(require.size() != 0) {
		printf("                require=%s\n",require.c_str());
		eval = new BLEvaluator();
	}

	BLManager::getObject()->registerSteppingAction(this);
	BLManager::getObject()->registerCallback(this,2);
}

void MaterialFilter::UserSteppingAction(const G4Step *step)
{
	// active in all states

	// only handle steps that enter a new volume
	G4VPhysicalVolume *preVol = 
				step->GetPreStepPoint()->GetPhysicalVolume();
	G4VPhysicalVolume *postVol = 
				step->GetPostStepPoint()->GetPhysicalVolume();
	if(preVol == postVol) return;

	// only handle this material (in the new volume)
	if(postVol->GetLogicalVolume()->GetMaterial() != material) return;

	G4Track *track = step->GetTrack();
	G4ParticleDefinition *def = track->GetDefinition();
	if(killSet.count(def) > 0 ||
			     (keepSet.size() > 0 && keepSet.count(def) == 0)) {
		track->SetTrackStatus(fStopAndKill);
		++nKilled;
		if(manager->getSteppingVerbose() > 0) 
			printf("material %s killed track\n",
				  postVol->GetLogicalVolume()->GetMaterial()->
				  			GetName().c_str());
	} else if(eval != 0) {
		eval->setTrackVariables(track,BLCOORD_GLOBAL);
		if(eval->evaluate(require) == 0.0) {
			track->SetTrackStatus(fStopAndKill);
			++nKilled;
			if(manager->getSteppingVerbose() > 0) 
				printf("material %s killed track\n",
				  postVol->GetLogicalVolume()->GetMaterial()->
				  			GetName().c_str());
		}
	}
}

void MaterialFilter::callback(int type)
{
	printf("material %s: %lu killed\n",material->GetName().c_str(),nKilled);
}

