//	g4blgui.cc

#include <stdio.h>

#include <QApplication>
#include <QString>
#include <QDir>

#include "G4blGuiWindow.h"
#include "../g4bl/Geant4Data.hh"
#include "../g4bl/Util.hh"

QString g4bl_dir;

int main(int argc, char *argv[])
{
	g4bl_dir = getG4blDir();
	QCoreApplication::addLibraryPath(g4bl_dir+"/lib");
	QCoreApplication::addLibraryPath(g4bl_dir+"/plugins");

	QApplication app(argc,argv);
	app.connect(&app,SIGNAL(lastWindowClosed()),&app,SLOT(quit()));

	QDir::setCurrent(QDir::homePath());

/*** no longer needed - using @rpath
#ifdef __linux__
	QString llp=g4bl_dir+"/lib:"+qgetenv("LD_LIBRARY_PATH").constData();
	qputenv("LD_LIBRARY_PATH",qPrintable(llp));
#endif
#ifdef __APPLE__
	QString llp=g4bl_dir+"/lib:"+qgetenv("DYLD_LIBRARY_PATH").constData();
	qputenv("DYLD_LIBRARY_PATH",qPrintable(llp));
#endif
***/
#ifdef WIN32
	// Need to add root's bin into PATH, so its DLLs are found.
	// (Root finds itself from where DLLs are loaded, not ROOTSYS.)
	// Fortunately, this program does not use Root.
	QString path=g4bl_dir + "\\root\\bin;" + qgetenv("PATH").constData();
	qputenv("PATH",qPrintable(path));
	printf("PATH=%s\n",qPrintable(path));
#endif

	if(!Geant4Data::setup()) Geant4Data::runG4bldata();

	new G4blGuiWindow();

	return app.exec();
}
