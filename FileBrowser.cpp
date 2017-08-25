#include "stdafx.h"
#include "FileBrowser.h"
#include <QDirModel>
#include <QHeaderView>
#include <QMessageBox>



class MyDirModel : public QDirModel
{
public:
	MyDirModel(QObject* parent) : QDirModel(parent) {}
	virtual int columnCount (const QModelIndex & parent = QModelIndex()) const
	{
	    if (parent.column() > 0)
		    return 0;
		return 1;
	}
};



FileBrowser::FileBrowser(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
	m_treemodel = new MyDirModel(this);
	m_treemodel->setFilter(QDir::AllDirs | QDir::NoDotAndDotDot);
	m_treemodel->setLazyChildCount(true);
	m_treemodel->setResolveSymlinks(false);
	m_treemodel->setReadOnly(true);

	m_listmodel = new QDirModel(this);
	m_listmodel->setFilter(QDir::Files);
	m_listmodel->setLazyChildCount(true);
	m_listmodel->setResolveSymlinks(false);
	m_listmodel->setReadOnly(true);

	QStringList names;
	names.append("*.off");
	names.append("*.obj");
	names.append("*.ply2");
	names.append("*.3ds");
	m_listmodel->setNameFilters(names);

	ui.dirView->setModel(m_treemodel);
	ui.dirView->header()->setVisible(false);

	QHeaderView *vh = ui.dirFiles->verticalHeader();
	vh->setVisible(false);
	vh->setResizeMode(QHeaderView::Fixed);
	vh->setDefaultSectionSize(17);

	ui.dirFiles->setModel(m_listmodel);
	ui.dirFiles->setColumnHidden(2, true);
	ui.dirFiles->setColumnHidden(3, true);
	ui.dirFiles->setSelectionBehavior(QAbstractItemView::SelectRows);

	QHeaderView *hh = ui.dirFiles->horizontalHeader();
	hh->setHighlightSections(false);
//	hh->resizeSection(0, ui.dirFiles->width() - 61);

	hh->resizeSection(1, 60);
	hh->setClickable(false);

	ui.dirFiles->setRootIndex(m_listmodel->index("C:\\"));


	connect(ui.dirView, SIGNAL(itemSelectionChanged()), this, SLOT(changedDir()));
	connect(ui.dirFiles, SIGNAL(doubleClicked(const QModelIndex&)), this, SLOT(selectFile(const QModelIndex&)));

	ProgSettings set;
	QString cdir = set.value("curDir", "").toString();
	if (!cdir.isEmpty())
		ui.dirView->setCurrentIndex(m_treemodel->index(cdir));
}

void FileBrowser::changedDir()
{
	const QModelIndex& index = ui.dirView->selectionModel()->currentIndex();
	if (!index.isValid())
		return;

	QString dir = m_treemodel->filePath(index);
	QModelIndex lind = m_listmodel->index(dir);
	ui.dirFiles->setRootIndex(lind);

	QHeaderView *hh = ui.dirFiles->horizontalHeader();
	hh->resizeSection(0, qMax(ui.dirFiles->width() - 81, 150));



	ProgSettings set;
	set.setValue("curDir", m_treemodel->filePath(ui.dirView->currentIndex()));
}

void FileBrowser::resizeEvent(QResizeEvent *event)
{
	QHeaderView *hh = ui.dirFiles->horizontalHeader();
	hh->resizeSection(0, qMax(ui.dirFiles->width() - 81, 150));
	QWidget::resizeEvent(event);
}

void FileBrowser::selectFile(const QModelIndex& index)
{
	QString f = m_listmodel->filePath(index);
	emit openFile(f);
	//QMessageBox::information(NULL, "", f);
}

