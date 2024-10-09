//	BLQt.h
/*
	Interface to the Geant4 Qt user interface and viewer.
	All references to Qt in the g4bl program are confined to BLQt.h
	and BLQt.cc.

	(This file is .h, not .hh, so the cmake AUTOMOC can find it; it has
	 no #ifdef G4BL_GUI because that would prevent the moc from working.
	 The #include-s of this file must be inside that #ifdef.)
*/

#ifndef BLQT_HH
#define BLQT_HH

#include <QString>
#include <QStringList>
#include <QToolBar>
#include <QLineEdit>
#include <QObject>

/// class BLQt is the interface to the Qt viewer and user interface.
class BLQt {
	static class BLQtToolBar *toolBar;
	static class G4UIExecutive *uie;
public:
	/// init1() does initialization for the Qt UI window, such as
	/// setting the plugin path and instantiating G4UIExecutive.
	static void init1();

	/// init2() re-organizes the Qt UI Window for G4beamline; it hides
	/// the Geant4 UI widgets and adds our toolbar.
	/// Once init2() has been called, call SessionStart() to
	/// open the Qt UI window and instantiate the viewer.
	static void init2();

	/// Session Start will create and open the Qt UI Window, and
	/// enter the Qt event loop.
	static void SessionStart();

	/// set the # events per image (= run)
	static void setEvPerImage(int evPerImage);

	/// returns the text for useful tips.
	static const char *usefulTips();
};

/// class BLQtToolBar is the G4beamline tool bar in the Qt user interface.
class BLQtToolBar : public QToolBar {
	Q_OBJECT
	QLineEdit *nEv;
	QStringList homeCommands;
	bool axesVisible;
	void getHomeCommands();
protected:
	void paintEvent(QPaintEvent *ev);
public:
	BLQtToolBar();
	void setNev(int n);
public slots:
	void run();
	void home();
	void setHome();
	void setViewpoint();
	void showAxes();
	void debug();
};

#endif // BLQT_HH
