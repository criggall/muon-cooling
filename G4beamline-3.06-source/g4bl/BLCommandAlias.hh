//	BLCommandAlias.hh

#ifndef BLCOMMANDALIAS_HH
#define BLCOMMANDALIAS_HH

#include "BLCommand.hh"

/**	class BLCommandAlias defines an alias for an existing command
 *
 *	Usage: To define "alias" as an aslias for "mycommand":
 *	class BLCMDmycommand : public BLCommand {
 *		...
 *	};
 *	BLCMDmycommand defaultMycommand;
 *	BLCommandAlias aliasMycommand("alias",defaultMycommand);
 *
 *	Note this is a global instantiation of BLCommandAlias, and that the
 *	second argument must be the global default instantiation of the desired
 *	command. This also works if BLCMDmycommand is derived from BLElement.
 *
 *	The alias will be registered as a command to BLCommand, and will get
 *	an appropriate entry in the help text.
 **/
class BLCommandAlias : public BLCommand {
	const char *alias;
	BLCommand *instance;
public:
	BLCommandAlias(const char *_alias, BLCommand &_instance) : BLCommand(){
		alias = strdup(_alias);
		instance = &_instance;
		registerCommand(instance->getCmdType());
		char tmp[256];
		sprintf(tmp,"Alias for '%s'.",instance->commandName().c_str());
		setSynopsis(strdup(tmp));
		setDescription("");
	}

	G4String commandName() { return G4String(alias); }

	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
		{ return instance->command(argv,namedArgs); }

};

#endif // BLCOMMANDALIAS_HH
