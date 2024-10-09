//	g4bldata.cc

#include <stdio.h>

#include <QApplication>
#include <QString>

#include "G4blDataWindow.h"
#include "../g4bl/Util.hh"

QString g4bl_dir;

int main(int argc, char *argv[])
{
	g4bl_dir = getG4blDir();
	QCoreApplication::addLibraryPath(g4bl_dir+"/lib");
	QCoreApplication::addLibraryPath(g4bl_dir+"/plugins");

	QApplication app(argc,argv);
	app.connect(&app,SIGNAL(lastWindowClosed()),&app,SLOT(quit()));

#ifdef __linux__
	QString llp=g4bl_dir + "/lib:" + qgetenv("LD_LIBRARY_PATH").constData();
	qputenv("LD_LIBRARY_PATH",qPrintable(llp));
#endif

	new G4blDataWindow();

	return app.exec();
}
