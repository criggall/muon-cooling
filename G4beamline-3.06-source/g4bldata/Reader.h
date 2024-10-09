//	Reader.h

#ifndef READER_H
#define READER_H

#include <QObject>
#include <QString>
#include <QByteArray>
#include <QFile>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>

/**	copyUrlToFile() will do just that.
 *
 *	Note that copyUrlToFile() "blocks" until the copy is complete or
 *	an error occurs. While "blocked", it loops calling 
 *	QCoreApplication::processEvents() and usleep(1000). Thus the user
 *	interface remains responsive, QWidget-s are updated, etc.
 *
 *	The URL should be complete; the file can be relative or absolute path.
 *	The url can be a "file:" URL, or can be a filename or path.
 *
 *	progress should point to a QObject with slots progress(qint64,qint64)
 *	and finished().
 *
 *	returns 0 on success, error code on error.
 *
 *	Prints error messages to stdout, but is silent on success.
 **/
int copyUrlToFile(QString url, QString file, QObject *progress=0);

/**	class Reader will read a file or URL.
 *
 *	Supports these URLs: http://..., ftp://..., file:/...
 *	Also filenames, either absolute or relative paths.
 *
 *	TYPICAL USAGE:
 *	Reader *r = Reader::create("/path/to/file"); // or "http://...", etc.
 *	connect(r,SIGNAL(finished()), ...);
 *	connect(r,SIGNAL(progress(qint64,qint64)), ...);
 *	for(;;) {
 *		QByteArray ba = r->read()
 *		if(ba.size() == 0 || r->error() != 0) break;
 *		//... use the data in ba
 *	}
 *	delete r;
 *
 *	Note that read() "blocks" until some data are available, the 
 *	end-of-file is reached, or an error occurs. While "blocked", it loops
 *	calling QCoreApplication::processEvents() and usleep(1000). Thus the 
 *	user interface remains responsive, QWidget-s are updated, etc.
 **/
class Reader : public QObject {
	Q_OBJECT
protected:
	int error_number;
	qint64 bytes;
	qint64 total;
public:
	Reader() : QObject(), error_number(0), bytes(0), total(0) { }
	virtual ~Reader() { }
public:
	/// Interprets fileOrUrl as a file or a URL, and constructs the
	/// corresponding FileReader or UrlReader.
	/// NOTE: the caller must delete the object.
	static Reader *create(QString fileOrUrl);

	/// Reads data as it becomes available. Returns an empty QByteArray
	/// at EOF or upon error. "Blocks" waiting for data, processing Qt
	/// events and sleeping for 1ms.
	virtual QByteArray read() = 0;

	/// returns zero if no error, otherwise returns the Qt error code.
	int error() { return error_number; }
signals:
	/// Emitted as data is transferred and read.
	void progress(qint64 bytes, qint64 total);

	/// Emitted when the data transfer is finished. No more data will come.
	void finished();
};

class FileReader : public Reader {
	Q_OBJECT
	QFile *file;
public:
	FileReader(QString filename);
	virtual ~FileReader() { if(file) file->close(); }
	virtual QByteArray read();
};

class UrlReader : public Reader {
	Q_OBJECT
	bool done;
	QNetworkAccessManager nam;
	QNetworkReply *reply;
	QByteArray downloadedData;
public:
	UrlReader(QString urlName);
	virtual ~UrlReader() { if(reply) delete reply; }
	virtual QByteArray read();
private slots:
	void dataReady();
	void downloadComplete();
};

#endif // READER_H
