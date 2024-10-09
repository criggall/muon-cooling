//	G4blGuiWindow.cc

#include <assert.h>

#include <QApplication>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QScrollBar>
#include <QTextCursor>
#include <QFileDialog>
#include <QFileInfo>
#include <QMenuBar>
#include <QAction>
#include <QUrl>
#include <QDesktopServices>
#include <QMessageBox>
#include <QRegExp>
#include <QList>

#include "G4blGuiWindow.h"
#include "../g4bl/Util.hh"
#include "../g4bl/Geant4Data.hh"

extern QString g4bl_dir;

G4blGuiWindow::G4blGuiWindow() : QMainWindow(0), process(0), abortTries(0),
				outputFile(0), outputStream(0), inViewer(false),
				recentFile(0)
{
	connect(QApplication::instance(),SIGNAL(aboutToQuit()), 
							this,SLOT(quit()));
	setWindowTitle("g4blgui");
	QWidget *w = new QWidget;
	setCentralWidget(w);
	QVBoxLayout *layout = new QVBoxLayout;
	w->setLayout(layout);
	layout->setSpacing(1);

	QMenu *fileMenu = menuBar()->addMenu("&File");
	QMenu *viewMenu = menuBar()->addMenu("&View");
	QMenu *toolMenu = menuBar()->addMenu("&Tools");
	QMenu *helpMenu = menuBar()->addMenu("&Help");
	fileMenu->addAction("Open File",this,SLOT(browseInputFile()));
	fileMenu->addSeparator();
	QMenu *rf = fileMenu->addMenu("Open Recent");
	fileMenu->addSeparator();
	fileMenu->addAction("Quit",this,SLOT(quit()));
	recentFile = new RecentFile(this,rf);
#ifdef G4BL_VISUAL
	viewMenu->addAction("Viewer = OpenGLStoredQt", this,SLOT(setVisual()));
	viewMenu->addAction("Viewer = Other Viewer", this,SLOT(setVisual()));
#endif
	viewMenu->addAction("No Viewer",this,SLOT(noVisual()));
	toolMenu->addAction("Do Regression Tests",
						this,SLOT(doRegressionTests()));
	helpMenu->addAction("Help",this,SLOT(setHelp()));
	helpMenu->addAction("Users Guide",this,SLOT(setHelp()));
	helpMenu->addAction("Validation Document",this,SLOT(setHelp()));
	helpMenu->addAction("README",this,SLOT(setHelp()));

	QHBoxLayout *lay = new QHBoxLayout;
	lay->addWidget(new QLabel("Input file:"));
	inputFile = new QLineEdit("");
	lay->addWidget(inputFile);
	QPushButton *b = new QPushButton("Browse");
	connect(b,SIGNAL(clicked()),this,SLOT(browseInputFile()));
	lay->addWidget(b);
	lay->addWidget(new QLabel("  Directory:"));
	directory = new QLineEdit("");
	lay->addWidget(directory);
	layout->addLayout(lay);

	lay = new QHBoxLayout;
	lay->addWidget(new QLabel("Parameters:"));
	parameters = new QLineEdit("");
	lay->addWidget(parameters);
	layout->addLayout(lay);

	lay = new QHBoxLayout;
#ifdef G4BL_VISUAL
	visual = new QCheckBox("Visualization (OpenGLStoredQt)");
	visual->setChecked(true);
	lay->addWidget(visual);
	evPerImageLabel = new QLabel("          Events per Image:");
	lay->addWidget(evPerImageLabel);
	evPerImage = new QSpinBox;
	evPerImage->setValue(1);
	connect(visual,SIGNAL(clicked()),this,SLOT(visClicked()));
	lay->addWidget(evPerImage);
#else
	visual = 0;
	evPerImageLabel = 0;
	evPerImage = 0;
	lay->addWidget(new QLabel("(Visualization not supported)"));
#endif
	lay->addStretch();
	helpOutput = new QPushButton("Show Help");
	connect(helpOutput,SIGNAL(clicked()),this,SLOT(doHelpOrOutput()));
	lay->addWidget(helpOutput);
	runAbort = new QPushButton("Run");
	connect(runAbort,SIGNAL(clicked()),this,SLOT(doRunOrAbort()));
	lay->addWidget(runAbort);
	layout->addLayout(lay);

	select = new QStackedWidget;
	helpBrowser = new QTextBrowser;
	select->addWidget(helpBrowser);
	outputText = new QPlainTextEdit;
	select->addWidget(outputText);
	layout->addWidget(select);
	outputText->setPlainText("The output from the most recent run will be "
		"displayed here.\n\nYou can switch between Help and Output at "
		"any time without affecting either.");
	QStringList helpPath; helpPath << g4bl_dir+"/share/g4beamline";
	helpBrowser->setSearchPaths(helpPath);
	helpBrowser->setSource(QUrl("file:Help.html"));
	connect(helpBrowser,SIGNAL(sourceChanged(const QUrl&)),
				this,SLOT(helpSourceChanged(const QUrl&)));

	doHelpOrOutput();
	resize(800,600);
	show();
}

QStringList G4blGuiWindow::getArgs()
{
	QStringList ret;

	ret << inputFile->text();
	QString line = parameters->text();
	while(line.length() > 0) {
		static QRegExp ws("^\\s*");
		static QRegExp nameq1("^([A-Za-z0-9_]+)=\"([^\"]*)\"");
		static QRegExp nameq2("^([A-Za-z0-9_]+)='([^']*)'");
		static QRegExp name1("^([A-Za-z0-9_]+)=([^\\s]*)\\s");
		static QRegExp name2("^([A-Za-z0-9_]+)=([^\\s]*)$");
		static QRegExp posq1("^\"([^\"]+)\"");
		static QRegExp posq2("^'([^']+)'");
		static QRegExp pos1("^([^\\s]+)\\s");
		static QRegExp pos2("^([^\\s]+)$");
		line = line.replace(ws,"");
		if(line.length() == 0) break;
		// remove quotes just like the shell, as g4bl expects
		if(line.contains(nameq1)) {
			line.replace(0,nameq1.matchedLength(),"");
			ret << nameq1.cap(1) + "=" + nameq1.cap(2);
		} else if(line.contains(nameq2)) {
			line.replace(0,nameq2.matchedLength(),"");
			ret << nameq2.cap(1) + "=" + nameq2.cap(2);
		} else if(line.contains(name1)) {
			line.replace(0,name1.matchedLength(),"");
			ret << name1.cap(1) + "=" + name1.cap(2);
		} else if(line.contains(name2)) {
			line.replace(0,name2.matchedLength(),"");
			ret << name2.cap(1) + "=" + name2.cap(2);
		} else if(line.contains(posq1)) {
			line.replace(0,posq1.matchedLength(),"");
			ret << posq1.cap(1);
		} else if(line.contains(posq2)) {
			line.replace(0,posq2.matchedLength(),"");
			ret << posq2.cap(1);
		} else if(line.contains(pos1)) {
			line.replace(0,pos1.matchedLength(),"");
			ret << pos1.cap(1);
		} else if(line.contains(pos2)) {
			line.replace(0,pos2.matchedLength(),"");
			ret << pos2.cap(1);
		} else {
			QMessageBox::warning(this,"G4beamline Error",
				"Syntax error in Parameters");
			break; // prevent infinite loop
		}
	}

#ifdef G4BL_VISUAL
	if(visual->isChecked()) {
		ret << "viewer=best," + QString::number(evPerImage->value());
	}
#endif

	return ret;
}

void G4blGuiWindow::openFile(QString filename)
{
	inputFile->setText(filename);
	QFileInfo fi(filename);
	// cd to the directory containing the input file
	QDir::setCurrent(fi.dir().absolutePath());
	recentFile->fileOpened(filename);
	outputText->clear();
}

void G4blGuiWindow::appendOutputText(QString _text)
{
	// remove benign but scary Root error messages (leave '\n')
	// (but message can be split if it comes first)
	static QString prepend;
	static bool first=true;
	QString text(prepend + _text);
	prepend = "";
	if(first && text.size() <  49) {
		prepend = text;
		return;
	}
	first = false;
	removeBenignErrorMessages(text);

	// moving the cursor and scrollbars will be invisible to the
	// user, as there is no opportunity to refresh the widget.
	// get current cursor and scrollbar positions
	int vpos = outputText->verticalScrollBar()->value();
	int vmax = outputText->verticalScrollBar()->maximum();
	int hpos = outputText->horizontalScrollBar()->value();
	QTextCursor cursor = outputText->textCursor();
	// insert new text at the end (append() adds a \n we don't want)
	outputText->moveCursor(QTextCursor::End);
	outputText->insertPlainText(text);
	// restore cursor
	if(cursor.atEnd())
		outputText->moveCursor(QTextCursor::End);
	else
		outputText->setTextCursor(cursor);
	// restore scroll bars, going to new max if was at old max
	if(vpos != vmax) {
		outputText->verticalScrollBar()->setValue(vpos);
		outputText->horizontalScrollBar()->setValue(hpos);
	} else {
		outputText->verticalScrollBar()->setValue(
		outputText->verticalScrollBar()->maximum());
		outputText->horizontalScrollBar()->setValue(
		outputText->horizontalScrollBar()->minimum());
	}
}

void G4blGuiWindow::removeBenignErrorMessages(QString &text)
{
	static const char *list[] = {
		// From Root 5
		"Error: cannot open file \"iostream\"[^\\n]*\\n",
		"\\*\\*\\* Interpreter error recovered[^\\n]*\\n",
		// From old Max OS X visualizations
		"QCocoaView handleTabletEvent:[^\\n]*\\n",
		// from Root 6
		"ERROR in cling::[^\\n]*\\n",
		"Invoking:\\n",
		"  LC_ALL=C /Library/Developer/CommandLineTools/[^\\n]*\\n",
		"Results was:\\n",
		"With exit code 256\\n",
		" *\\^[^\\n]*\\n",
		"^input_line_[^\\n]* fatal error:[^\\n]*\\n",
		"^#include[^\\n]*\\n",
		"^In file included from[^\\n]*\\n",
		"^C:\\Program Files[^\\n]*\\n",
		"^ *struct _CrtEnableIf[^\\n]*\\n",
		"^ *typedef struct[^\\n]*\\n",
		"^ *included multiple times[^\\n]*\\n",
		"^Warning in <TClassTable[^\\n]*\\n",
	};
	static QList<QRegExp> regexp;
	if(regexp.isEmpty()) {
		for(unsigned i=0; i<sizeof(list)/sizeof(list[0]); ++i) {
			regexp.append(QRegExp(list[i]));
		}
	}

	for(int i=0; i<regexp.size(); ++i)
		text.remove(regexp[i]);
}

void G4blGuiWindow::quit()
{
	if(process) process->kill();
	::exit(99);
}

void G4blGuiWindow::browseInputFile()
{
	QString filename = QFileDialog::getOpenFileName(this,
				"Open G4beamline input",
				QDir::currentPath(),
				"G4bl input (*.g4bl *.in);;All files (*)");
	if(!filename.isEmpty()) {
		openFile(filename);
	}
}

void G4blGuiWindow::doRunOrAbort()
{
	QString text = runAbort->text();
	if(text == "Run") {
		if(process != 0) {
		    QMessageBox::warning(this,"G4beamline Gui",
			"A process is running, let it complete, or Abort it.");
		    return;
		}
		abortTries = 0;
		inViewer = false;
		outputText->clear();
		select->setCurrentIndex(1);
		helpOutput->setText("Show Help");
		QString dir = directory->text();
		if(!dir.isEmpty()) {
			if(dir.contains("/") || dir.contains("\\")) {
				appendOutputText("*** Invalid directory -- "
					"Must be a simple directory name.");
				return;
			}
			// cd to the directory containing the input file
			QFileInfo fi(inputFile->text());
			QDir::setCurrent(fi.dir().absolutePath());
			// create or use the directory
			QDir d=QDir::current();
			d.mkdir(dir) ; // let it fail if it already exists
			if(!d.cd(dir) || 
			   !QDir::setCurrent(d.absolutePath())) {
				appendOutputText("*** Cannot create directory '"
					+dir+"'");
				return;
			}
		}
		outputFile = new QFile("g4bl.out");
		if(!outputFile->open(QIODevice::WriteOnly)) {
			appendOutputText("*** Cannot create 'g4bl.out'\n\n");
		}
		outputStream = new QTextStream(outputFile);
		runAbort->setText("Abort");
		process = new QProcess(this);
		process->setProcessChannelMode(QProcess::MergedChannels);
		connect(process,SIGNAL(readyReadStandardOutput()),
			this,SLOT(updateProcessOutput()));
		connect(process,SIGNAL(finished(int,QProcess::ExitStatus)),
			this,SLOT(processFinished(int,QProcess::ExitStatus)));
		connect(process,SIGNAL(error(QProcess::ProcessError)),
			this,SLOT(error(QProcess::ProcessError)));
		QString path = g4bl_dir + "/bin/g4bl";
		QStringList args = getArgs();
		process->start(path,args);
		if(!dir.isEmpty()) QDir::current().cdUp();
	} else if(text == "Abort") {
		if(!process) return;
		// (avoid extraneous output)
		disconnect(process,SIGNAL(error(QProcess::ProcessError)),0,0);
		if(inViewer) {
			appendOutputText("\n--KILLED while in Viewer\n");
			process->kill();
			inViewer = false;
		} else if(abortTries++ == 0) {
			appendOutputText("\n--ABORTED (will exit at end of "
				"event, NTuples are short but valid)\n");
			process->terminate();
			// let output drain and finished() be called
		} else {
			// first one didn't take -- kill it
			appendOutputText("\n--KILLED (will exit immediately, "
				"Ntuples and Root files are invalid)\n");
			process->kill();
		}
	} else if(text == "Abort Tests") {
		if(!process) return;
		// (avoid extraneous output)
		disconnect(process,SIGNAL(error(QProcess::ProcessError)),0,0);
		if(abortTries++ == 0) {
			appendOutputText("\n--ABORTED\n");
			process->terminate();
			// let output drain and finished() be called
		} else {
			// first one didn't take -- kill it
			appendOutputText("\n--KILLED\n");
			process->kill();
		}
	}
}

void G4blGuiWindow::doHelpOrOutput()
{
	QString text = helpOutput->text();
	if(text == "Show Help") {
		select->setCurrentIndex(0);
		helpOutput->setText("Show Output");
	} else {
		select->setCurrentIndex(1);
		helpOutput->setText("Show Help");
	}
}

void G4blGuiWindow::visClicked()
{
#ifdef G4BL_VISUAL
	if(visual->isChecked()) {
		evPerImageLabel->setVisible(true);
		evPerImage->setVisible(true);
	} else {
		evPerImageLabel->setVisible(false);
		evPerImage->setVisible(false);
	}
#endif
}

void G4blGuiWindow::helpSourceChanged(const QUrl &src)
{
	if(src.path() == "/G4blData") {
		Geant4Data::runG4bldata();
		helpBrowser->home();
	}
}

void G4blGuiWindow::setVisual()
{
#ifdef G4BL_VISUAL
	QAction *a = dynamic_cast<QAction*>(sender());
	if(a != 0) {
	    QString s = a->text();
	    if(s.contains("OpenGL") || s.contains("OGL")) {
		visual->setChecked(true);
		visClicked();
	    } else {
		visual->setChecked(false);
		visClicked();
	    	QMessageBox::information(this,"Other Viewer",
			"The OpenGLStoredQt viewer is generally the best "
			"one to use, as it provides a full-featured 3-D viewer "
			"controlled by the mouse.\n\n"
			"To use some other viewer, add a parameter "
			"'viewer=name,N,M' to the parameters. Here name is the "
			"name of the viewer, N is the number of events per "
			"image, and M is the number of images to generate "
			"(N and M default to 1 if omitted).\n\n"
			"The following viewers are usually available:\n"
			"  ASCIITree (ATree)\n"
			"  DAWNFILE (DAWNFILE)\n"
			"  G4HepRep (HepRepXML)\n"
			"  G4HepRepFile (HepRepFile)\n"
			"  OpenGLImmediateX (OGLIX)\n"
			"  OpenGLImmediateXm (OGLI, OGLIXm)\n"
			"  OpenGLStoredX (OGLSX)\n"
			"  OpenGLStoredQt (OGLSQt)\n"
			"  OpenGLImmediateQt (OGLIQt)\n"
			"  RayTracer (RayTracer)\n"
			"  RayTracerX (RayTracerX)\n"
		);
	    }
	}
#endif
}

void G4blGuiWindow::noVisual()
{
#ifdef G4BL_VISUAL
	visual->setChecked(false);
	visClicked();
#endif
}

void G4blGuiWindow::setHelp()
{
	QAction *a = dynamic_cast<QAction*>(sender());
	QUrl url;
	url.setScheme("file");
	if(a != 0) {
	    QString s = a->text();
	    if(s == "Help") {
		select->setCurrentIndex(0);
		helpOutput->setText("Show Output");
	    } else if(s == "Users Guide") {
		url.setPath(g4bl_dir+"/doc/G4beamlineUsersGuide.pdf");
		if(!QDesktopServices::openUrl(url))
			QMessageBox::warning(this,"G4beamline Gui",
				"Cannot open URL "+url.toString());
	    } else if(s == "Validation Document") {
		url.setPath(g4bl_dir+"/doc/G4beamlineValidation.pdf");
		if(!QDesktopServices::openUrl(url))
			QMessageBox::warning(this,"G4beamline Gui",
				"Cannot open URL "+url.toString());
	    } else if(s == "README") {
		url.setPath(g4bl_dir+"/doc/README.txt");
		if(!QDesktopServices::openUrl(url))
			QMessageBox::warning(this,"G4beamline Gui",
				"Cannot open URL "+url.toString());
	    }
	}
}

void G4blGuiWindow::doRegressionTests()
{
	if(process != 0) {
		QMessageBox::warning(this,"G4beamline Gui",
			"A process is running, let it complete, or Abort it.");
		return;
	}

	outputText->clear();
	outputFile = 0;
	outputStream = 0;
	abortTries = 0;

	select->setCurrentIndex(1);
	helpOutput->setText("Show Help");
	runAbort->setText("Abort Tests");

	process = new QProcess(this);
	process->setProcessChannelMode(QProcess::MergedChannels);
	connect(process,SIGNAL(readyReadStandardOutput()),
		this,SLOT(updateProcessOutput()));
	connect(process,SIGNAL(finished(int,QProcess::ExitStatus)),
		this,SLOT(processFinished(int,QProcess::ExitStatus)));
	connect(process,SIGNAL(error(QProcess::ProcessError)),
		this,SLOT(error(QProcess::ProcessError)));
	QString path = g4bl_dir + "/bin/g4bltest";
	QStringList args;
	args << "--loop";
	process->start(path,args);
}

void G4blGuiWindow::error(QProcess::ProcessError error)
{
	switch(error) {
	case QProcess::FailedToStart:
		appendOutputText("-- CANNOT START\n");
		break;
	case QProcess::Crashed:
		appendOutputText("-- CRASHED\n");
		break;
	default:
		appendOutputText("-- UKNOWN ERROR\n");
		break;
	}
	// processFinished() will be called very soon
}

void G4blGuiWindow::updateProcessOutput()
{
	assert(process != 0);
	QByteArray data = process->readAllStandardOutput();
	QString s(data);
	if(s.isEmpty()) return;
	int in = s.lastIndexOf("You have entered a viewer secondary X");
	int out = s.lastIndexOf("Secondary X event loop exited.");
	if(in >= 0 && in > out) inViewer = true;
	if(out >= 0 && out > in) inViewer = false;
	appendOutputText(s);
	if(outputStream != 0) *outputStream << s;
}

void G4blGuiWindow::processFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
	assert(process != 0);
	updateProcessOutput();
	QString s("\nProgram Exited normally.\n");
	if(exitCode != 0)
		s = QString("\n--ERROR: Program exited with error code %1\n").
							arg(exitCode);
	if(exitStatus == QProcess::CrashExit)
		s = QString("\n--ERROR: Program crashed\n");
	appendOutputText(s);
	process->deleteLater();
	process = 0;
	if(outputFile != 0) outputFile->close();
	if(outputStream != 0) delete outputStream;
	outputStream = 0;
	if(outputFile != 0) delete outputFile;
	outputFile = 0;
	runAbort->setText("Run");
}

