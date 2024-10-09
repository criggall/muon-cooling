//	BLQt.cc

#ifdef G4BL_VISUAL

#ifdef G4BL_GUI

#ifndef G4UI_USE_QT
#define G4UI_USE_QT
#endif

#include <QTextEdit>
#include <QIcon>
#include <QApplication>
#include <QLabel>
#include <QFormLayout>
#include <QWidget>
#include <QSizePolicy>
#include <QAction>
#include <QMainWindow>
#include <QGLWidget>
#include <QList>
#include <QLayout>
#include <QDockWidget>
#include <QToolBar>
#include <QTimer>
#include <QCoreApplication>

#include "G4UImanager.hh"
#include "G4VViewer.hh"
#include "G4Scene.hh"
#include "G4ViewParameters.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4UIQt.hh"

#include "BLQt.h"
#include "BLParam.hh"
#include "BLVisManager.hh"

extern G4String g4bl_dir;   // in g4bl.cc
extern G4String input_file; // in g4bl.cc

BLQtToolBar * BLQt::toolBar = 0;
G4UIExecutive * BLQt::uie = 0;

void printQtChildren(const QObject *w, const char *prefix="") 
{
	printf("%s%s %s\n",prefix,w->metaObject()->className(),
						qPrintable(w->objectName()));
	QString s(prefix);
	s += "        ";
	Q_FOREACH(const QObject *child,w->children()) {
		printQtChildren(child,qPrintable(s));
	}
}

void BLQt::init1()
{
	// add our lib and plugins to Qt's library search path
	QString dir(g4bl_dir);
	QCoreApplication::addLibraryPath(dir+"/lib");
	QCoreApplication::addLibraryPath(dir+"/plugins");

	// save current UI Session and create a G4UIQt session window (which
	// automatically becomes the current UI session)
	G4UIsession *orgSession = G4UImanager::GetUIpointer()->GetSession();
	static char *argv[] = {(char *)"g4bl",0};
	uie = new G4UIExecutive(1,argv,"Qt");

	// send G4cout and G4cerr to stdout and stderr (not G4UIQt)
	G4UImanager::GetUIpointer()->SetCoutDestination(orgSession);
}

void BLQt::init2()
{
	// reorganize the main window to hide all the Geant4 command stuff
	// (unused in G4beamline), and add our toolbar

	G4UIQt *session = dynamic_cast<G4UIQt*>(
				G4UImanager::GetUIpointer()->GetSession());
	BLAssert(session != 0);
	QMainWindow *main = session->GetMainWindow();
	BLAssert(main != 0);

	// Hide the Geant4 QDockWidget-s (not needed in G4beamline)
	QList<QDockWidget*> docks = main->findChildren<QDockWidget*>("");
	Q_FOREACH(QDockWidget *d, docks) {
		d->hide();
	}

	// hide the Geant4 QToolBar-s (not needed in G4beamline)
	// NOTE: This is too early, so BLQtToolbar does it in paintEvent()

	// add the G4beamline toolbar on the left, and resize to a better size
	toolBar = new BLQtToolBar();
	main->addToolBar(Qt::LeftToolBarArea, toolBar);
	QString dir(g4bl_dir);
       QApplication::setWindowIcon(QIcon(dir+"/share/g4beamline/g4blicon.png"));
	main->setWindowTitle("G4beamline   "+QString(input_file));
	main->resize(700,500);

	// change the text in the "useful tips" tab to reflect the changes.
	QList<QTextEdit*> edits = main->findChildren<QTextEdit*>("");
	edits[0]->setHtml(usefulTips());

	// Display the Qt Widget hierarchy, if parameter "QtWidgets" is defined
	if(Param.isDefined("QtWidgets")) {
		printf("Qt Widget hierarchy, from Main Window:\n");
		printQtChildren(main);
	}
	fflush(stdout);
}

void BLQt::SessionStart()
{
	if(uie) uie->SessionStart();
}

void BLQt::setEvPerImage(int evPerImage)
{
	if(toolBar) toolBar->setNev(evPerImage);
}

BLQtToolBar::BLQtToolBar() : QToolBar(), axesVisible(false)
{
	nEv = new QLineEdit("1");
	nEv->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
	nEv->setToolTip("Number of events per run.");
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

	addSeparator();
	a = addAction("Debug");
	connect(a, SIGNAL(triggered()), this, SLOT(debug()));
	a->setToolTip("Perform some debugging action (see code).");

	QTimer::singleShot(500,this,SLOT(setHome()));
}

void BLQtToolBar::setNev(int n)
{
	nEv->setText(QString::number(n));
}

void BLQtToolBar::run()
{
	int n = nEv->text().toInt();
	BLRunManager::getObject()->BeamOn(n);
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
	if(!axesVisible)
	    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/add/axes");
	axesVisible = true;
}

void BLQtToolBar::debug()
{
	printf("BLQtToolBar::debug -- does nothing\n"); fflush(stdout);
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

void BLQtToolBar::paintEvent(QPaintEvent *ev)
{
	static bool first=true;
	if(first) {
		// hide the Geant4 QToolBar-s (not needed in G4beamline)
		G4UIQt *session = dynamic_cast<G4UIQt*>(
				G4UImanager::GetUIpointer()->GetSession());
		BLAssert(session != 0);
		QMainWindow *main = session->GetMainWindow();
		QList<QToolBar*> toolbars = main->findChildren<QToolBar*>("");
		Q_FOREACH(QToolBar *t, toolbars) {
			if(t == this) continue;
			t->hide();
		}
		first = false;
	}

	QToolBar::paintEvent(ev);
}

const char *BLQt::usefulTips()
{
// nested macros are needed to stringize the value of a preprocessor symbol
#define QUOTE(A) QUOTE_(A)
#define QUOTE_(A) #A

	return "&nbsp;<BR/>\n"
	"<H2><center>G4beamline " QUOTE(G4BLVERSION) "</center></H2>\n"
	"<P>This is the OpenGLStoredQt viewer. It displays a 3-D view of the "
	"world in the main portion of the window, with a tool bar on the left. "
	"The number at the top is the number of events per run (image).</P>\n"
	"<P>In the 3-D viewer, the following actions are available:<BR/>\n"
	"<table><tr><td>Right Click</td><td>Context Menu</td></tr>\n"
	"<tr><td>Middle Click</td><td>(No Op)</td></tr>\n"
	"<tr><td>Left Click and Drag</td><td>Rotate or Pan (See Context Menu, default is Rotate)</td></tr>\n"
	"<tr><td>(Shift) Click and Drag</td><td>Pan</td></tr>\n"
	"<tr><td>Middle Click</td><td>(No Op)</td></tr>\n"
	"<tr><td>Arrow Keys (U,D,L,R)</td><td>Pan (move U,D,L,R)</td></tr>\n"
	"<tr><td>(Shift) Arrow Keys</td><td>Rotate</td></tr>\n"
	"<tr><td>(Alt) Arrow Keys</td><td>Rotate</td></tr>\n"
	"<tr><td>(Cmd) + -</td><td>Zoom In and Out</td></tr>\n"
	"<tr><td>(Alt) + -</td><td>Increase/Decrease rotation for arrow keys</td></tr>\n"
	"</table>\n"
	"<P>When the context menu says the Mouse Action is Rotate (the "
	"default), then (Shift) Click and Drag will pan the display.</P>\n"
	"<P>Picking of volumes or events is not yet implemented.</P>\n"
	"<P>At the top of the toolbar is a small region where a right-click "
	"will bring up a menu that lets you select 'Scene tree, Help, History' "
	"or 'Output' -- these are advanced toolbars that assist you in using "
	"Geant4 commands; un-checking the bottom line will hide the toolbar, "
	"and you cannot get it back.</P>\n"
	;
}

#endif // G4BL_GUI

#else  // !G4BL_VISUAL
int dummy_g4bltoolbar=0;
#endif // G4BL_VISUAL
