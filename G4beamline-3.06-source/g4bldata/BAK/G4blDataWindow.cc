//	G4blDataWindow.cc

const char *instructions =
 "This program will download Geant4 datasets from MuonsInc and install them.\n"
 "Normally the BaseURL and Destination should be left unchanged.\n"
 "\n"
 "The physics lists used in your simulations determine which data sets are required.\n"
 "If you're not sure which you need, and have the space, download them all."
 "\n";

const char *changeDestinationInstructions =
 "G4beamline looks for the datasets in $HOME/Geant4Data. On windows it also\n"
 "looks in C:\\Geant4Data. If there is a file '.data' in the installation "
 "directory,\nit contains the absolute path of the Geant4Data directory.\n"
 "\n"
 "If you use a different location, you must create the '.data' file.\n";

#include <assert.h>

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QPushButton>
#include <QDialog>
#include <QDir>
#include <QCoreApplication>
#include <QProcess>
#include <QMessageBox>

#include "G4blDataWindow.h"
#include "Reader.h"
#include "../g4bl/Geant4Data.hh"

#ifdef WIN32
#define EXT ".exe"
#else
#define EXT ""
#endif

extern QString g4bl_dir;

G4blDataWindow::G4blDataWindow() : QMainWindow(), checkBoxes(),
			progressLabels(), currentProgressLabel(0),
			currentLength(), first(true)
{
	connect(QCoreApplication::instance(),SIGNAL(aboutToQuit()), 
							this,SLOT(abort()));
	setWindowTitle("g4bldata");
	QWidget *w = new QWidget;
	setCentralWidget(w);
	QVBoxLayout *layout = new QVBoxLayout;
	w->setLayout(layout);
	layout->setSpacing(1);

	QLabel *label=new QLabel("g4bldata - Download datasets for G4beamline");
	label->setStyleSheet("font-weight:bold; font-size:20pt; "
							"text-align:center;");
	layout->addWidget(label);
	notice = new QLabel("");
	layout->addWidget(notice);
	label = new QLabel(instructions);
	layout->addWidget(label);

	QHBoxLayout *lay = new QHBoxLayout;
	lay->addWidget(new QLabel("    Base URL:"));
	baseURL = new QLineEdit(Geant4Data::baseURL.c_str());
	lay->addWidget(baseURL);
	layout->addLayout(lay);

	lay = new QHBoxLayout;
	lay->addWidget(new QLabel("Destination:"));
	destination = new QLabel(qgetenv("HOME")+"/Geant4Data");
	lay->addWidget(destination);
	lay->addStretch();
	QPushButton *b = new QPushButton("Change");
	connect(b,SIGNAL(clicked()),this,SLOT(changeDestination()));
	lay->addWidget(b);
	layout->addLayout(lay);

	QFormLayout *form = new QFormLayout;
	form->setVerticalSpacing(2);
#ifndef STUB
	std::vector<Geant4Data> list = Geant4Data::list();
	form->addRow(new QLabel(""));
	label = new QLabel("Required for all physics lists:");
	label->setStyleSheet("font-weight:bold;");
	form->addRow(label);
	for(std::vector<Geant4Data>::iterator i=list.begin(); i!=list.end(); ++i) {
		if(!i->required) continue;
		directories << i->dir.c_str();
		lengths << QString::number(i->dl_size);
		checkBoxes << new QCheckBox(i->name.c_str());
		progressLabels << new QLabel("???");
		checkBoxes.last()->setFixedWidth(300);
		form->addRow(checkBoxes.last(),progressLabels.last());
	}
	form->addRow(new QLabel(""));
	label = new QLabel("Additional Datasets:");
	label->setStyleSheet("font-weight:bold;");
	form->addRow(label);
	for(std::vector<Geant4Data>::iterator i=list.begin(); i!=list.end(); ++i) {
		if(i->required) continue;
		directories << i->dir.c_str();
		lengths << "0";
		checkBoxes << new QCheckBox(i->name.c_str());
		progressLabels << new QLabel("???");
		checkBoxes.last()->setFixedWidth(300);
		form->addRow(checkBoxes.last(),progressLabels.last());
	}
#else // STUB
	for(const char **p=fileList; (*p)!=0; ++p) {
		QString s(*p);
		if(s.endsWith(":")) {
// less whitespace	form->addRow(new QLabel(""));
			label = new QLabel(s);
			label->setStyleSheet("font-weight:bold;");
			form->addRow(label);
		} else if(s.endsWith("/")) {
			directories << s;
		} else if(s[0].isDigit()) {
			lengths << s;
		} else {
			checkBoxes << new QCheckBox(s);
			progressLabels << new QLabel("???");
			checkBoxes.last()->setFixedWidth(300);
			form->addRow(checkBoxes.last(),progressLabels.last());
		}
	}
#endif // STUB
	assert(checkBoxes.size() == progressLabels.size());
	assert(checkBoxes.size() == directories.size());
	assert(checkBoxes.size() == lengths.size());
	layout->addLayout(form);
// less whitespace	layout->addWidget(new QLabel(""));

	lay = new QHBoxLayout;
	b = new QPushButton("Download and Install Checked Datasets");
	b->setDefault(true);
	connect(b,SIGNAL(clicked()),this,SLOT(go()));
	lay->addWidget(b);
	cancel = new QPushButton("Cancel");
	connect(cancel,SIGNAL(clicked()),this,SLOT(abort()));
	lay->addWidget(cancel);
	layout->addLayout(lay);

	changeDestination(); // first time no dialog, just scan destination

	resize(100,100);
	show();
	raise();
}

int G4blDataWindow::extractAll(QString filename)
{
	// find the tar program
	QString tar(g4bl_dir+"/bin/bsdtar" EXT);
	if(!QFile::exists(tar)) tar = g4bl_dir+"/bin/tar" EXT;
	if(!QFile::exists(tar)) tar = "/bin/tar" EXT;
	if(!QFile::exists(tar)) tar = "/usr/bin/tar" EXT;
	// run tar to extract the files
	QProcess p(this);
	p.setProcessChannelMode(QProcess::ForwardedChannels);
	QStringList args;
	args << "-xzf" << filename;
	p.start(tar,args);
	int timeLimit = 300000; // milliseconds
	if(!p.waitForFinished(timeLimit) || p.exitCode() != 0 ||
	   			p.exitStatus() != QProcess::NormalExit) {
		return -1;
	}
	return 0;
}

void G4blDataWindow::go()
{
	QString url = baseURL->text();
	if(!url.endsWith("/")) url += "/";
	QString dest = destination->text();
	if(!dest.endsWith("/")) dest += "/";
	if(!QDir(dest).exists())
		QDir(dest).mkpath(".");
	QString current = QDir::currentPath();
	if(!QDir::setCurrent(dest)) {
		printf("Cannot cd to '%s'\n",qPrintable(dest));
		QMessageBox::critical(this,"Error","Cannot cd to '"+dest+"'");
		return;
	}
	for(int i=0; i<checkBoxes.size(); ++i) {
		if(!checkBoxes[i]->isChecked()) continue;
		checkBoxes[i]->setChecked(false);
		currentProgressLabel = progressLabels[i];
		currentLength = lengths[i];
		currentProgressLabel->setText("Starting...");
		QString file = checkBoxes[i]->text() + ".tar.gz";
		fprintf(stderr,"downloading %s\n",qPrintable(url+file));
		if(copyUrlToFile(url+file,dest+file,this) != 0) {
			currentProgressLabel->setText("ERROR DOWNLOADING");
			fprintf(stderr,"ERROR DOWNLOADING\n");
			continue;
		}
		fprintf(stderr,"Extracting...\n");
		currentProgressLabel->setText("Extracting...");
		currentProgressLabel->repaint();
		if(extractAll(dest+file) != 0) {
			currentProgressLabel->setText("ERROR EXTRACTING");
			continue;
		}
		fprintf(stderr,"Complete\n");
		currentProgressLabel->setText("Complete");
 		QFile::remove(dest+file);
	}
	QDir::setCurrent(current);
	currentProgressLabel = 0;
	cancel->setText("Close");
}

void G4blDataWindow::progress(qint64 bytes, qint64 total)
{
	if(!currentProgressLabel) return;
	QString v = QString::number(bytes) + "/" + QString::number(total);
	if(total < 0) v = QString::number(bytes) + "/" + currentLength;
	currentProgressLabel->setText(v);
}

void G4blDataWindow::finished()
{
	//printf("finished\n");
}

void G4blDataWindow::changeDestination()
{
	if(!first) {
		QDialog *d = new QDialog(this);
		d->setWindowTitle("g4bldata - Change Destination");
		QVBoxLayout *layout = new QVBoxLayout;
		d->setLayout(layout);
		layout->addWidget(new QLabel(changeDestinationInstructions));
		QHBoxLayout *lay = new QHBoxLayout;
		lay->addWidget(new QLabel("New Destination:"));
		QLineEdit *newDest = new QLineEdit(destination->text());
		lay->addWidget(newDest);
		layout->addLayout(lay);
		lay = new QHBoxLayout;
		lay->addStretch();
		QPushButton *b = new QPushButton("Save");
		connect(b,SIGNAL(clicked()),d,SLOT(accept()));
		lay->addWidget(b);
		b = new QPushButton("Cancel");
		connect(b,SIGNAL(clicked()),d,SLOT(reject()));
		lay->addWidget(b);
		lay->addStretch();
		layout->addLayout(lay);
		d->show();
		if(d->exec() == QDialog::Accepted) {
			destination->setText(newDest->text());
		}
		delete d;
	}
	first = false;

	// scan destination for existing dataset directories
	int n=0;
	QDir dir(destination->text());
	for(int i=0; i<checkBoxes.size(); ++i) {
		if(dir.exists(directories[i])) {
			progressLabels[i]->setText("Present");
			checkBoxes[i]->setChecked(false);
			++n;
		} else {
			progressLabels[i]->setText("");
			bool v = !checkBoxes[i]->text().startsWith("G4TENDL");
			checkBoxes[i]->setChecked(v);
		}
	}

	if(n >= 6) {
		notice->setText("\nThe Geant4Data directory has been found.\n"
			"You may not need to download anything.\n");
		notice->setStyleSheet("font-weight:bold; color:red;");
	}
}

void G4blDataWindow::abort()
{
	::exit(99);
}
