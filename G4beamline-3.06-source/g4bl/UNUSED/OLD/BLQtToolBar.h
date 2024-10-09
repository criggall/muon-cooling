//	BLQtToolBar.hh

#ifndef BLQTTOOLBAR_HH
#define BLQTTOOLBAR_HH

#include <QString>
#include <QStringList>
#include <QToolBar>
#include <QLabel>
#include <QLineEdit>
#include <QFormLayout>
#include <QWidget>
#include <QSizePolicy>

class BLQtToolBar : public QToolBar {
	Q_OBJECT
	QLineEdit *nEv;
	QStringList homeCommands;
	bool axesVisible;
	void getHomeCommands();
public:
	BLQtToolBar();
public slots:
	void run();
	void home();
	void setHome();
	void setViewpoint();
	void showAxes();
};

#endif // BLQTTOOLBAR_HH
