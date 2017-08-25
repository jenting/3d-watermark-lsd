#pragma once

#include <qobject.h>
#include <qevent.h>
#include <qfile.h>
#include <qtextstream.h>
#include <qdialog.h>
#include <q3dockwindow.h>
#include <q3textedit.h>

class StatusDialog : public Q3DockWindow/*QDialog*/ {
public:
	StatusDialog(QWidget* parent = 0, const char* name = 0) :
	  Q3DockWindow/*QDialog*/(Q3DockWindow::InDock, parent, name) {		
	}

protected:
	void resizeEvent(QResizeEvent* resizeEvent);
};

class QGlobalLogEvent : public QEvent {
private:
	QStringList m_strings;
	int m_logLevel;
public:
	QGlobalLogEvent(const int logLevel, const QStringList& strings) : 
		QEvent(QEvent::User), m_strings(strings), m_logLevel(logLevel) {

	}

	QGlobalLogEvent(const int logLevel, const QString& string) : 
		QEvent(QEvent::User), m_logLevel(logLevel) {
		m_strings.append(string);
	}

	const QStringList& strings() const {
		return m_strings;
	}

	const int logLevel() const {
		return m_logLevel;
	}
};

#define SDFLOG(logLevel, s) GlobalLog::log(logLevel, s)
#define SDFLOG10(s) GlobalLog::log(10,s)
#define SDFLOG8(s) GlobalLog::log(8,s)
#define SDFLOG6(s) GlobalLog::log(6,s)
#define SDFLOG4(s) GlobalLog::log(4,s)
#define SDFLOG2(s) GlobalLog::log(2,s)

class GlobalLog : public QObject {
	Q_OBJECT
protected:
	static GlobalLog* m_log;
	static const int RECENT_STRINGS_SIZE = 50;

protected:
	QFile m_file;
	QTextStream* m_stream;
	int m_logLevel;	
	//QStringList m_recent;

	StatusDialog* m_statusDialog;
	//QListBox* m_listbox;
	Q3TextEdit* m_logBox;

	/** construct a log object, protected so nonoe can construct it */
	GlobalLog(QWidget* parent, const QString& logFile, const int logLevel);

	~GlobalLog();

	void writeToFile(const QStringList& strings);

	void writeToScreen(const QStringList& strings);

	void createStatusDialog();

	void updateStatusDialog();
	
public:
	/** create a global log object */
	static void init(QWidget* parent, const QString& logFile, const int logLevel);

	/** release global log object and close everything */
	static void destroy();

	static void logm(const int logLevel, const QStringList& strings);
	static void log(const int logLevel, const QString& string);

	static const int logLevel() {
		return getLog()->m_logLevel;
	}

	/** get pointer to global log */
	static GlobalLog* getLog() {
		return m_log;
	}	

	/** override event function when log postings arrive */
	bool event(QEvent* e);

public slots:
	/** show log window */
	void ShowLog();
	/** hide log window */
	void HideLog();
};