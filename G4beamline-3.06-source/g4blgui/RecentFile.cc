//	RecentFile.cc

#include <QAction>

#include "RecentFile.h"
#include "G4blGuiWindow.h"

RecentFile::RecentFile(G4blGuiWindow *_window, QMenu *_menu) :
			window(_window), menu(_menu), recent(), settings(0)
{
	settings = new QSettings("muonsinc","g4blgui");
	recent = settings->value("RecentFiles").toStringList();
	fileOpened(QString());
}

void RecentFile::fileOpened(QString filename)
{
	if(!filename.isEmpty()) {
		if(!recent.isEmpty() && filename == recent.first())
			return;
		int j = recent.indexOf(filename);
		if(j >= 0) recent.removeAt(j);
		recent.prepend(filename);
		while(recent.size() > MAXRECENT)
			recent.removeLast();
		settings->setValue("RecentFiles",recent);
		settings->sync();
	}
	menu->clear();
	for(int i=0; i<recent.size(); ++i) {
		menu->addAction(recent[i],this,SLOT(recentFile()));
	}
}

void RecentFile::recentFile()
{
	QAction *s = dynamic_cast<QAction*>(sender());
	if(s == 0) return;
	QString f = s->text();
	if(f.isEmpty()) return;
	window->openFile(f); // calls fileOpened()
}
