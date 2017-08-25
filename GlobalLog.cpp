#include "stdafx.h"

#include "GlobalLog.h"

#include <qfileinfo.h>
#include <qdatetime.h>
#include <qapplication.h>
#include <qlayout.h>
#include <qdir.h>
#include <qfont.h>
#include <qmainwindow.h>
//Added by qt3to4:
#include <QResizeEvent>
#include <QEvent>

GlobalLog* GlobalLog::m_log = NULL;

void StatusDialog::resizeEvent(QResizeEvent* resizeEvent)
{
	ProgSettings ps;
	ps.beginGroup("status dialog");
	ps.writeEntry("width", width());
	ps.writeEntry("height", height());
	ps.writeEntry("x", x());
	ps.writeEntry("y", y());
	ps.endGroup();
}

void GlobalLog::init(QWidget* parent, const QString& logFile, const int logLevel)
{
	m_log = new GlobalLog(parent, logFile, logLevel);
}

void GlobalLog::destroy()
{
	delete m_log;
}

void GlobalLog::logm(const int logLevel, const QStringList& strings)
{
	if (logLevel>=getLog()->m_logLevel) {
		QGlobalLogEvent* gle = new QGlobalLogEvent(logLevel, strings);
		QApplication::postEvent(getLog(), gle); 
	}
}

void GlobalLog::log(const int logLevel, const QString& string)
{
	if (logLevel>=getLog()->m_logLevel) {
		QGlobalLogEvent* gle = new QGlobalLogEvent(logLevel, string);
		QApplication::postEvent(getLog(), gle);
	}	
}

GlobalLog::GlobalLog(QWidget* parent, const QString& logFile, const int logLevel):
	QObject(parent),
	m_file(logFile), 
	m_logLevel(logLevel), 
	m_stream(NULL), 
	m_statusDialog(NULL)
{
	QFileInfo fi(m_file);
	if (fi.exists() && fi.size()>1024*1024) {
		QDir dir(fi.filePath());
		dir.rename(m_file.name(), m_file.name()+".old");
		dir.remove(logFile);
	}

	if (m_file.open(QIODevice::WriteOnly | QIODevice::Append)) {
		m_stream = new QTextStream(&m_file);

	}

	if (!m_statusDialog) createStatusDialog();
}

GlobalLog::~GlobalLog()
{

	if (m_file.isOpen()) {
		m_file.close();
	}

	if (m_stream) 
		delete m_stream;

}

void GlobalLog::writeToFile(const QStringList& strings)
{
	QString dt = QDateTime::currentDateTime().toString("[yyyy-MM-dd hh:mm:ss.zzz] ");

	QStringList::const_iterator it = strings.begin();
	for (;it != strings.end();it++) {
		*m_stream << dt << *it << "\r\n";
	}
}

void GlobalLog::writeToScreen(const QStringList& strings)
{
	for (QStringList::const_iterator it = strings.begin(); it != strings.end(); it++) {
		m_logBox->append(*it);
	}
}

bool GlobalLog::event(QEvent* e)
{
	if (e->type() != QEvent::User) return false;

	QGlobalLogEvent* gle = (QGlobalLogEvent*)e;

	writeToScreen(gle->strings());

	if (m_stream && gle->logLevel() >= m_logLevel)
		writeToFile(gle->strings());

	if (m_statusDialog && m_statusDialog->isVisible())
		updateStatusDialog();

	return true;
}

/*void GlobalLog::createStatusDialog()
{
	if (m_statusDialog) return;

	ProgSettings ps;
	ps.beginGroup("status dialog");
	int width = ps.readNumEntry("width", 120);
	int height = ps.readNumEntry("height", 600);
	int x = ps.readNumEntry("x", 0);
	int y = ps.readNumEntry("y", 0);
	ps.endGroup();

	m_statusDialog = new StatusDialog((QWidget*)parent());
	m_statusDialog->setCaption("SDF Log");
	m_statusDialog->resize(width, height);
	if (x) m_statusDialog->move(x,y);

	QVBoxLayout* vbox = new QVBoxLayout(m_statusDialog);

	//m_listbox = new QListBox(m_statusDialog);
	m_logBox = new QTextEdit(m_statusDialog);
	m_logBox->setTextFormat(Qt::LogText);
	m_logBox->setMaxLogLines(RECENT_STRINGS_SIZE);
	QStyleSheetItem* header = new QStyleSheetItem(m_logBox->styleSheet(), "h");
	header->setColor("blue");
	header->setFontWeight(QFont::Bold);
	header->setFontUnderline(true);	

	vbox->addWidget(m_logBox);	
}*/

void GlobalLog::createStatusDialog()
{
	if (m_statusDialog) return;

	ProgSettings ps;
	ps.beginGroup("status dialog");
	int width = ps.readNumEntry("width", 120);
	int height = ps.readNumEntry("height", 600);
	int x = ps.readNumEntry("x", 0);
	int y = ps.readNumEntry("y", 0);
	ps.endGroup();

	m_statusDialog = new StatusDialog((QWidget*)parent());
	m_statusDialog->setVerticallyStretchable(true);
	m_statusDialog->setResizeEnabled(true);
	
//	((QMainWindow*)parent())->moveDockWindow(m_statusDialog, Qt::DockRight);
		
	m_logBox = new Q3TextEdit(m_statusDialog);
	m_logBox->setTextFormat(Qt::LogText);
	m_logBox->setMaxLogLines(RECENT_STRINGS_SIZE);
	Q3StyleSheetItem* header = new Q3StyleSheetItem(m_logBox->styleSheet(), "h");
	header->setColor("blue");
	header->setFontWeight(QFont::Bold);
	header->setFontUnderline(true);	

	m_statusDialog->boxLayout()->addWidget(m_logBox);	
	m_logBox->show();

	m_statusDialog->hide();
}

void GlobalLog::updateStatusDialog()
{
	if (!m_statusDialog || !m_logBox) return;

	//m_logBox->clear();
	//for (QStringList::iterator it = m_recent.begin(); it != m_recent.end(); it++)
		//m_logBox->append(*it);	
	//m_logBox->setCurrentItem(m_listbox->count());
	//m_logBox->ensureCurrentVisible();
}

void GlobalLog::ShowLog()
{
	if (!m_statusDialog) createStatusDialog();

	if (!m_statusDialog) return;

	if (m_statusDialog->isVisible()) {
		m_statusDialog->hide();
	} else {
		//updateStatusDialog();
		m_statusDialog->show();
	}	
}

void GlobalLog::HideLog()
{
	if (!m_statusDialog) return;

	m_statusDialog->hide();

}