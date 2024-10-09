//	G4blGuiWindow.h

#ifndef G4BLGUIWINDOW_H
#define G4BLGUIWINDOW_H

#include <stdio.h>

#include <QMainWindow>
#include <QString>
#include <QStringList>
#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>
#include <QSpinBox>
#include <QStackedWidget>
#include <QPlainTextEdit>
#include <QTextBrowser>
#include <QProcess>
#include <QFile>
#include <QTextStream>
#include <QMenu>

#include "RecentFile.h"

class G4blGuiWindow : public QMainWindow {
	Q_OBJECT
	QLineEdit *inputFile;
	QLineEdit *directory;
	QLineEdit *parameters;
	QCheckBox *visual;
	QLabel *evPerImageLabel;
	QSpinBox *evPerImage;
	QPushButton *helpOutput;
	QPushButton *runAbort;
	QStackedWidget *select;
	QPlainTextEdit *outputText;
	QTextBrowser *helpBrowser;
	QProcess *process;
	int abortTries;
	QFile *outputFile;
	QTextStream *outputStream;
	bool inViewer;
	RecentFile *recentFile;
public:
	G4blGuiWindow();
	QStringList getArgs();
	void openFile(QString filename);
	void appendOutputText(QString text);
	void removeBenignErrorMessages(QString &text);
public slots:
	void quit();
	void browseInputFile();
	void doRunOrAbort();
	void doHelpOrOutput();
	void visClicked();
	void helpSourceChanged(const QUrl &src);
	void setVisual();
	void noVisual();
	void setHelp();
	void doRegressionTests();
private slots:
	void error(QProcess::ProcessError error);
	void updateProcessOutput();
	void processFinished(int exitCode, QProcess::ExitStatus exitStatus);
};

#endif // G4BLGUIWINDOW_H
