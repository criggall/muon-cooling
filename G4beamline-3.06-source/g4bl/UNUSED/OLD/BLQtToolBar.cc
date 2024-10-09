//	BLQtToolBar.cc

#ifdef G4BL_VISUAL

#ifdef G4BL_GUI

#include <QAction>

#include "G4UImanager.hh"
#include "G4VViewer.hh"
#include "G4Scene.hh"
#include "G4ViewParameters.hh"

#include "BLQtToolBar.h"
#include "BLVisManager.hh"

BLQtToolBar::BLQtToolBar() : QToolBar(), axesVisible(false)
{
	nEv = new QLineEdit("1");
	nEv->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
	nEv->setToolTip("Number of events per run.");
	//nEv->setMinimumWidth(40);
	//QSizePolicy sp = nEv->sizePolicy();
	//sp.setHorizontalStretch(1);
	//nEv->setSizePolicy(sp);
	addWidget(nEv);
	QAction *a = addAction("Start Run");
	connect(a, SIGNAL(triggered()), this, SLOT(run()));
	a->setToolTip("Start a run with the above number of events; draw them.");

	addSeparator();
	a = addAction("Home");
	connect(a, SIGNAL(triggered()), this, SLOT(home()));
	a->setToolTip("Restore the viewer to Home position.");
	a = addAction("Set Home");
	connect(a, SIGNAL(triggered()), this, SLOT(setHome()));
	a->setToolTip("Set the Home position to the current view.");

	addSeparator();
	a = addAction("Top View");
	connect(a, SIGNAL(triggered()), this, SLOT(setViewpoint()));
	a->setToolTip("Show system viewed from top (+y).");
	a = addAction("Left View");
	connect(a, SIGNAL(triggered()), this, SLOT(setViewpoint()));
	a->setToolTip("Show system viewed from left side (+x).");
	a = addAction("Right View");
	connect(a, SIGNAL(triggered()), this, SLOT(setViewpoint()));
	a->setToolTip("Show system viewed from right side (-x).");
	a = addAction("Front View");
	connect(a, SIGNAL(triggered()), this, SLOT(setViewpoint()));
	a->setToolTip("Show system viewed from front (+z).");
	a = addAction("Rear View");
	connect(a, SIGNAL(triggered()), this, SLOT(setViewpoint()));
	a->setToolTip("Show system viewed from rear (-z).");
	a = addAction("Bottom View");
	connect(a, SIGNAL(triggered()), this, SLOT(setViewpoint()));
	a->setToolTip("Show system viewed from bottom (-y).");

	addSeparator();
	a = addAction("Show Axes");
	connect(a, SIGNAL(triggered()), this, SLOT(showAxes()));
	a->setToolTip("Display (Red,Green,Blue) Axes at the origin.");

	getHomeCommands();
}

void BLQtToolBar::run()
{
	int n = nEv->text().toInt();
	QString cmd("/run/beamOn "+QString::number(n));
	G4UImanager::GetUIpointer()->ApplyCommand(qPrintable(cmd));
}

void BLQtToolBar::home()
{
	Q_FOREACH(QString s, homeCommands) {
		G4UImanager::GetUIpointer()->ApplyCommand(qPrintable(s));
	}
}

void BLQtToolBar::setHome()
{
	getHomeCommands();
}

void BLQtToolBar::setViewpoint()
{
	QAction *a = dynamic_cast<QAction*>(sender());
	if(a != 0) {
		QString text = a->text();
		QStringList cmds;
		// make it unlikely new viewpoint is along up vector
		cmds << "/vis/viewer/set/upVector 0.315 -.823 .51";
		if(text.startsWith("Front"))
			cmds << "/vis/viewer/set/viewpointVector 0 0 1"
			    << "/vis/viewer/set/upVector 0 1 0";
		else if(text.startsWith("Rear"))
			cmds << "/vis/viewer/set/viewpointVector 0 0 -1"
			    << "/vis/viewer/set/upVector 0 1 0";
		else if(text.startsWith("Top"))
			cmds << "/vis/viewer/set/viewpointVector 0 1 0"
			    << "/vis/viewer/set/upVector 1 0 0";
		else if(text.startsWith("Bottom"))
			cmds << "/vis/viewer/set/viewpointVector 0 -1 0"
			    << "/vis/viewer/set/upVector -1 0 0";
		else if(text.startsWith("Left"))
			cmds << "/vis/viewer/set/viewpointVector 1 0 0"
			    << "/vis/viewer/set/upVector 0 1 0";
		else if(text.startsWith("Right"))
			cmds << "/vis/viewer/set/viewpointVector -1 0 0"
			    << "/vis/viewer/set/upVector 0 1 0";
		Q_FOREACH(QString s, cmds) {
		    G4UImanager::GetUIpointer()->ApplyCommand(qPrintable(s));
		}
	}
}

void BLQtToolBar::showAxes()
{
	G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/add/axes");
}

void BLQtToolBar::getHomeCommands()
{
	G4Point3D p = BLVisManager::getObject()->GetCurrentScene()
						->GetStandardTargetPoint();
	G4VViewer *viewer = BLVisManager::getObject()->GetCurrentViewer();
	const G4ViewParameters &param = viewer->GetViewParameters();
	G4String s = param.CameraAndLightingCommands(p).c_str();
	homeCommands = QString(s.c_str()).split("\n");
}

#endif // G4BL_GUI

#else  // !G4BL_VISUAL
int dummy_g4bltoolbar=0;
#endif // G4BL_VISUAL
