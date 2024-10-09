//	G4blDataWindow.h

#ifndef G4BLDATAWINDOW_H
#define G4BLDATAWINDOW_H

#include <stdio.h>

#include <QMainWindow>
#include <QLineEdit>
#include <QLabel>
#include <QCheckBox>
#include <QList>
#include <QPushButton>

class G4blDataWindow : public QMainWindow {
	Q_OBJECT
	QLineEdit *baseURL;
	QLabel *destination;
	QLabel *notice;
	QList<QCheckBox*> checkBoxes;
	QList<QLabel*> progressLabels;
	QList<QString> directories;
	QList<QString> lengths;
	QLabel *currentProgressLabel;
	QString currentLength;
	QPushButton *cancel;
	bool first;
public:
	G4blDataWindow();
	virtual ~G4blDataWindow() { } // QWidgets are deleted automatically
	int extractAll(QString filename);
	void closeEvent(QCloseEvent *ev) { abort(); }
public slots:
	void go();
	void progress(qint64 bytes, qint64 total);
	void finished();
	void changeDestination();
	void abort();
};

#endif // G4BLDATAWINDOW_H
