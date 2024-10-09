//	Reader.cc

#include <QCoreApplication>

#include "Reader.h"

int copyUrlToFile(QString url, QString file, QObject *progress)
{
	int retval = -99;
	FILE *out = fopen(qPrintable(file),"wb");
	if(!out) {
		printf("Cannot Write '%s'\n",qPrintable(file));
		return -1;
	}
	Reader *reader = Reader::create(url);
	if(progress != 0) {
		progress->connect(reader,SIGNAL(progress(qint64,qint64)),
				progress,SLOT(progress(qint64,qint64)));
		progress->connect(reader,SIGNAL(finished()),
				progress,SLOT(finished()));
	}
	qint64 bytesRead=0;
	for(;;) {
		QByteArray ba = reader->read();
		retval = reader->error();
		if(ba.size() == 0 || retval != 0) break;
		bytesRead += ba.size();
		size_t n = fwrite(ba.data(),1,ba.size(),out);
		if(n != ba.size()) {
			printf("Write error on '%s'\n",qPrintable(file));
			delete reader;
			fclose(out);
			return -1;
		}
	}
	fclose(out);
	delete reader;
	return retval;
}

Reader *Reader::create(QString fileOrUrl)
{
	if(fileOrUrl.startsWith("http:") || fileOrUrl.startsWith("ftp:") ||
	   fileOrUrl.startsWith("file:"))
		return new UrlReader(fileOrUrl);
	return new FileReader(fileOrUrl);
}

FileReader::FileReader(QString filename) : Reader()
{
	file = new QFile(filename);
	if(!file->open(QIODevice::ReadOnly)) {
		error_number = file->error();
		delete file;
		file = 0;
	}
	if(file) total = file->size();
}

QByteArray FileReader::read()
{
	QCoreApplication::processEvents();
	QByteArray ret;
	if(file) {
		ret = file->read(4000);
		error_number = file->error();
		bytes += ret.size();
		emit progress(bytes,total);
	}
	if(ret.size() == 0) emit finished();
	return ret;
}

UrlReader::UrlReader(QString urlName) : Reader(), done(false), nam(),
						reply(0), downloadedData()
{
	QUrl url(urlName);
	QNetworkRequest request(url);
	reply = nam.get(request);
	reply->setTextModeEnabled(false);
	connect(reply,SIGNAL(readyRead()),this,SLOT(dataReady()));
	connect(reply,SIGNAL(error(QNetworkReply::NetworkError)),
					this,SLOT(dataReady()));
	connect(reply,SIGNAL(finished()),this,SLOT(downloadComplete()));
	connect(reply,SIGNAL(downloadProgress(qint64,qint64)),
					this,SIGNAL(progress(qint64,qint64)));
	dataReady(); // in case of glare (data ready before connect completes)
}

QByteArray UrlReader::read()
{
	QByteArray ret;
	for(;;) {
		QCoreApplication::processEvents();
		if(downloadedData.size() > 0 || error_number != 0 || done) break;
		//usleep(1000);
	}
	ret = downloadedData.left(4000);
	downloadedData.remove(0,ret.size());
	if(ret.size() == 0)
		emit finished();
	return ret;
}

void UrlReader::dataReady()
{
	QByteArray ba = reply->readAll();
	downloadedData.append(ba);
	error_number = reply->error();
}

void UrlReader::downloadComplete()
{
	done = true;
}

